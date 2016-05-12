import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
from operator import itemgetter

from Bio import SeqIO
from collections import Counter
import csv
import itertools
import glob, os


import random_colors

PFAM_IDS = ["AHH", "Colicin-DNase", "DUF1524", "DUF968", "Endonuclea_NS_2", "Endonuclease_1", 
	    "Endonuclease_7", "Endonuclease_NS", "GH-E", "HNH", "HNH_2", "HNH_3", "HNH_4",
	    "HNH_5", "HNHc_6", "ICEA", "LHH", "MH1", "NinG", "RecA_dep_nuc", "WHH", "zf-His_Me_endon", 
	    "Tox-HNH-EHHH","Tox-HNH-HHH","DUF1364","Tox-SHH", "RE_Alw26IDE"]

PFAM_IDS_2 = ["PF01844","PF01223","PF03165","PF13391","PF14279","PF07510","PF04231","PF13392","PF13395",
	      "PF13930","PF14412","PF02945","PF14414","PF12639","PF05766","PF16784","PF14411",
	      "PF06147","PF05551","PF14410","PF16786","PF05315", "PF09665"]
 
gis_dir = "gis"
WORKING_DIR = "."
all_seqs_file = "all_seqs_5_70_nolimits.fas"
pdb_mappings = "pdb_mapping.txt" # niefiltrowane przez GRDB
whole_seqs = "hits_round_2_70_ws.fas"
ws_pfam_mapping = "hits_round_2_70_ws_pfam.txt"
pfam_mapping = "hits_round_2_70_pfam_mapping_wholegis.txt"
cog_mapping = "all_seqs_5_70_cog_mapping.txt"
kog_mapping = "all_seqs_5_70_kog_mapping.txt"
gis_with_limits = os.path.join(WORKING_DIR, "all_seqs_5_70_gis_limits.txt")
final_cluster_dir = "/home/jagoda/Desktop/annotation_test"
limits = "hits_round_2_70_gis_and_limits.txt"
cd_hits = "hits_round_2_70.fas.clstr"
CDD_DB = "/home/jagoda/DB/cddid_all.tbl"
EASY_GO = os.path.join(WORKING_DIR, "gene_ontology.txt")
taxonomy_file = "taxonomy_round_2.txt"
PROTEIN_NAMES = "protein_names.txt"
FUNCTIONS = "domains_list.txt"


def get_gis_dict(gis_dir):
	"""
	Make gis dict in the following format cl_name:[gis_list]
	"""
	gis_dict = {}
	clusters = glob.glob(os.path.join(gis_dir,"*.txt"))
	for cl in clusters:
		cl_name = os.path.basename(cl).split(".")[0]
		gis_list = open(cl,"r").read().splitlines()
		stripped = []
		for gi in gis_list:
			stripped.append(gi.rstrip())
		gis_dict[cl_name] = stripped
	
	return gis_dict

def gather_sequences(fasta, gis_dict):
	"""
	gather all sequences belonging to clusters
	"""
	dir_path = os.path.join(WORKING_DIR, "sequences")

	record_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

	if not os.path.exists(dir_path):
    		os.makedirs(dir_path)
	
	for cluster in gis_dict:
		seq_file = open(os.path.join(dir_path, cluster+".fas"),"w")
		for gi in gis_dict[cluster]:
			if gi:							
				seq_file.write(">"+gi+"\n")
				seq_file.write(str(record_dict[gi].seq)+"\n")
			
		seq_file.close()
											
def get_taxonomy_count(taxonomy_file, gis_dict):
	"""
	Get taxonomy counts for every cluster.
	GI   taxa_name
	Draw a piechart
	"""
	taxonomies = open(taxonomy_file, "r").read().splitlines()
	dir_path = os.path.join(WORKING_DIR, "taxonomy")
	
	pie_colors = ['yellowgreen', 'mediumpurple', 'lightskyblue', 'lightcoral']

	if not os.path.exists(dir_path):
    		os.makedirs(dir_path)

	for cluster in gis_dict:
		tax_list = []
		for tax in taxonomies:
			if tax.split()[0] in gis_dict[cluster]:				
				tax_list.append(tax.split()[1])

		counted = dict(Counter(tax_list))
		
		num_taxa = sum(counted.values())
		sizes = [(float(x)/float(num_taxa))*100 for x in counted.values()]

		plt.pie(sizes,   
        		colors=pie_colors,
        		shadow=False,
        		startangle=70)

		plt.axis('equal')

		plt.savefig(cluster+"_tax.png", dpi=300)
		plt.close()
	
def get_pdb90(gi_dict, pdb_mappings):
	"""
	Map clusters to pdb files
	"""

	pdbs = open(pdb_mappings, "r").readlines()
	dir_path = os.path.join(WORKING_DIR, "pdb")
	
	cl2pdb = {}

	if not os.path.exists(dir_path):
    		os.makedirs(dir_path)
	
	for cluster in gi_dict:
		cl2pdb[cluster] = []
		pdb_list_file = open(os.path.join(WORKING_DIR,"pdb",cluster+".txt"), "w")
		for pb in pdbs:
			gi = pb.split()[1]
			if gi in gi_dict[cluster]:
				cl2pdb[cluster].append(pb.split()[0])
				pdb_list_file.write(pb.split()[0]+"\n")
	
		pdb_list_file.close()
	
	return cl2pdb

def map2limits(gis_limits):
	gis2limits = {}
	limits = open(gis_limits, "r").read().splitlines()
	
	for limit in limits:
		try:
			gis2limits[limit.split()[0]] = (int(limit.split()[1].split("-")[0]), int(limit.split()[1].split("-")[1]))
		except: pass
	return gis2limits

def domain_architectures(ws_pfam_mapping, gi_dict):
	"""
	Determine domain architectures for every seq.
	"""
	pfam_results = get_gi2pfams(ws_pfam_mapping)
	dir_path = os.path.join(WORKING_DIR, "architectures")
	
	if not os.path.exists(dir_path):
    		os.makedirs(dir_path)
	
	clusters_domains = {}
	for cluster in gi_dict:
		clusters_domains[cluster] = []
		for gi in gi_dict[cluster]:
			if gi in pfam_results:
				out = pfam_results[gi]
				clusters_domains[cluster].append(out)
	
	#return clusters_domains
	return filter_archs(clusters_domains)#,clusters_domains
	
def get_all_archs(unfiltered_dict):
	all_archs = {}
	all_archs["all"] = []
	for cluster in unfiltered_dict:
		for arch in unfiltered_dict[cluster]:
			all_archs["all"].append(arch)
	
	return all_archs	

def filter_archs(clusters_domains):
	"""
	Remove redundancy in architectures - repeats, multiple types, 1 aa domains
	"""
	
	final_archs = {}
	final_final_archs = {}	

	for cluster in clusters_domains:
		repeated = {}
		final_archs[cluster] = []
		for architectures in clusters_domains[cluster]:
			if tuple(set([x[2] for x in architectures])) not in repeated:
				final_archs[cluster].append(architectures)
			if tuple(set([x[2] for x in architectures])) not in repeated:
				repeated[tuple(set([x[2] for x in architectures]))] = 0
			repeated[tuple(set([x[2] for x in architectures]))] += 1
		
		for arch in final_archs[cluster]:
			key = tuple(set([x[2] for x in arch]))
			try:
				arch.append(repeated[key])	
			except:
				arch.append(0)

			if cluster not in final_final_archs:
				final_final_archs[cluster] = []

			final_final_archs[cluster].append(arch)

	final_final_archs = remove_repeats(final_final_archs)

	return final_final_archs
	
def remove_repeats(clusters_domains):
	"""
	Remove repeated domains, eg. RHS_repeat and repeat it with single box.
	"""
	new_clusters = {}

	for cluster in clusters_domains:	
		new_clusters[cluster] = []
		for architectures in clusters_domains[cluster]:
			new_archs = []
			added = []
			for domain in architectures:
				try:
					if domain[3] not in added:
						new_archs.append(domain)
						added.append(domain[3])
				except: 
						count = domain
						new_archs.append(count)

			new_clusters[cluster].append(new_archs)
	
	#return new_clusters			
	return sort_by_counts(new_clusters)

def sort_by_counts(clusters_dict):
	"""
	Sort architectures list by occurrence
	"""
	new_clusters = {}

	for cluster in clusters_dict:
			new_clusters[cluster] = []
			new_clusters[cluster] = sorted(clusters_dict[cluster], key=itemgetter(-1),reverse=False)
	
	return new_clusters

def overlap( (s1, e1), (s2, e2) ):
        return max(0, min(e1, e2) - max(s1, s2) + 1)

def get_gi2pfams(file, evalue_thresh=1e-5, max_overlap=0.5):
	"""
	Whole sequences pfam scan results.
	"""
        gi2pfams = {}
	gis2limits = map2limits(limits)
	gis = []
	
        for line in open(file):
                if line.startswith("#"):
                        continue
                tab = line.split()
                if len(tab) < 19:
                        continue
                gi = tab[3]
                gis.append(gi)
                
                if '|' in gi:
                        gi = gi.split('|')[1]
                pfam = tab[0]
		pfam_id = tab[1].split(".")[0]
                evalue = float(tab[6])

                if evalue > evalue_thresh:
                        continue
		
		if pfam in PFAM_IDS:
			continue

		if gi not in gi2pfams:
                        gi2pfams[gi] = []


                start, end = int(tab[17]), int(tab[18])
                length = end - start + 1
		
                for (start1, end1, _, _) in gi2pfams[gi]:
                        length1 = end1 - start1 + 1
                        n = min(length, length1)
                        if 1. * overlap((start1, end1), (start, end)) / n > max_overlap:
                                break
               	else:
                	gi2pfams[gi].append((start, end, pfam, pfam_id))


        for gi in gi2pfams:
		try:
			gi2pfams[gi].append((gis2limits[gi][0], gis2limits[gi][1], "His-Me", "His-Me"))
                except: pass
                gi2pfams[gi] = sorted(gi2pfams[gi])
	
	
	count = 0
	for gin in gi2pfams:
		if gin in gi2pfams:
			count+=1
				
        return gi2pfams

def run_mafft():
	"""
	Runs mafft on every cluster.
	"""
	dir_path = os.path.join(final_cluster_dir, "sequences")
	seq_files = glob.glob(os.path.join(dir_path,"*.fas"))

	for seq_file in seq_files:
		fname = os.path.basename(seq_file).split(".")[0]+".aln"
		cmd = "/opt/Mafft/bin/mafft-linsi --thread 60 %s > %s" % (seq_file, fname)
		os.system(cmd)

def map2colors(gi_dict, n_clusters=51):
	"""
	Returns color for each cluster.
	"""
	colors =  random_colors.generate_N_pastels(n_clusters)
	outfile = open(os.path.join(WORKING_DIR, "clusters_colors.txt"), "w")
	
	for num,cluster in enumerate(gi_dict):
		for gi in gi_dict[cluster]:
			outfile.write(gi+"\t"+colors[num]+"\n")

	outfile.close()	

def draw_architectures(arch_dict):
	"""
	Draw domain architectures using matplotlib.
	"""
	z = 0.3 #fixed width
	c = 0.1 #fixed height

	colors = functions2colors()

	for cluster in arch_dict:
		num_domains = len(arch_dict[cluster])
		max_len = max([len(x) for x in arch_dict[cluster]])
		fig6 = plt.figure()
		ax6 = fig6.add_subplot(111, aspect='equal')
		y = 0
		id = 100000

		lefts = []
		rights = []

		for arch in arch_dict[cluster]:	
			for i,dom in enumerate(arch[:-1]): #bad bad jagoda :(
				if dom[2] == "His-Me":
					id = i
					lefts.append(float(id))
					rights.append(float(len(arch[id:])))

		middle = (min(lefts)+max(rights))/2+0.5
		left = min(lefts)
		right = max(rights)
		# Get ID of His-Me domain
		for arch in arch_dict[cluster]:
			
			for i,dom in enumerate(arch[:-1]): #bad bad jagoda :(
				if dom[2] == "His-Me":
					id = i

			print "left",left
			print "right",right
			
			#middle = (left+right)/2+0.5
			print "middle", middle
			ax6.add_patch(patches.Rectangle((middle,y),z,c,facecolor="lightgrey",linewidth=0))
			
			ax6.text(0.5*(middle+middle+z), 0.5*(y+y+c), 'His-Me',
        			horizontalalignment='center', verticalalignment='center',
        			fontsize=4, color='black')
			x = middle-0.5
			for domain in arch[:id]:
				ax6.add_patch(patches.Rectangle((x,y),z,c,facecolor=colors[domain[3]],linewidth=0))
				
				ax6.text(0.5*(x+x+z), 0.5*(y+y+c), domain[2],
        			horizontalalignment='center', verticalalignment='center',
        			fontsize=4, color='black')
        			
				x-=0.3

			x = middle+0.5
			for domain in arch[id+1:-1]:
				ax6.add_patch(patches.Rectangle((x,y),z,c,facecolor=colors[domain[3]],linewidth=0))
				
				ax6.text(0.5*(x+x+z), 0.5*(y+y+c), domain[2],
        			horizontalalignment='center', verticalalignment='center',
        			fontsize=4, color='black')
				
				x+=0.3			

			# dodaj zliczenie
			ax6.text(0.49*(x+x+z), 0.5*(y+y+c), arch[-1],
        			horizontalalignment='center', verticalalignment='center',
        			fontsize=5, color='black')
			y+=0.1
			

		plt.ylim([0,num_domains*0.1])
		plt.xlim([left,right])
		plt.axis('off')
		#plt.subplots_adjust(left=0.001, right=0.002, top=0.002, bottom=0.001)
		#plt.savefig("domains_poster.pdf", format='pdf')

		fig6.savefig(cluster+'.png', dpi=150, bbox_inches='tight')	

	#pp.close()
	
def domains2file(arch_dict):
	"""
	Save domains to file.
	"""
	for cluster in arch_dict:
		domains_list = []
		outfile2 = open(cluster+"_domains_summary.txt", "w")
		outfile = open(cluster+"_domains.csv","w")
		for domain in arch_dict[cluster]:
			pfams = [x[2] for x in domain[:-1]]
			domains_list+=pfams
			outfile.write(",".join([str(x).strip('()').replace(","," ").replace("'","") for x in domain])+"\n")
		
		counter = Counter(domains_list).most_common()
		for cnt in counter:
			outfile2.write(cnt[0]+"\t"+str(cnt[1])+"\n")
			
		outfile2.close()
		outfile.close()
	
def domain2stats(unfiltered_domains):
	
	dom_counts = {}
	for cluster in unfiltered_domains:
		for arch in unfiltered_domains[cluster]:
			for domain in arch[:-1]:
				if domain[2] not in dom_counts:
					dom_counts[domain[2]] = 0
				dom_counts[domain[2]] += 1
	
def get_all_sequences(clusters_file, gis_dict):
	"""
	Get gis of ALL sequences (not 70%) 
	"""
	cd_hits = open(clusters_file, "r").readlines()

	cl_dict = {}
	lst = []
	cl_rep = "blabla"
	gis= []
	
	
	all_gis = {}

	for cluster in gis_dict:
		
		all_cl_gis = []
		gis = gis_dict[cluster]

		for line in cd_hits:
			if line.startswith(">Cluster"):
				if cl_rep in gis:
				        all_cl_gis+=lst
				lst = []
				cl_rep = ""

			else:
				lst.append(line.split()[2][:-3].split("_")[0][1:])
				if line.split()[3] == "*":
				        cl_rep = line.split()[2][:-3].split("_")[0][1:]

		all_gis[cluster] = all_cl_gis

	return all_gis

def get_pfam_mapping(gis_dict, pfam_mapping_file):
	"""
	Get PFAM ids for each cluster.
	"""
	pfams = open(pfam_mapping_file, "r").read().splitlines()
	pfams_dict = {}	
	cl2pfam = {}	

	for pf in pfams:
		pfams_dict[pf.split()[0]] = pf.split()[1]
	
	for cluster in gi_dict:
		cl2pfam[cluster] = []
		for gi in gi_dict[cluster]:
			if gi:	
				if gi in pfams_dict:
					cl2pfam[cluster].append(pfams_dict[gi])

		nonred = list(set(cl2pfam[cluster]))
		cl2pfam[cluster] = nonred

	return cl2pfam

def make_table(pdbs, pfams, cogs, kogs):
	"""
	Save table with data to csv.
	"""
	cluster_rows = sorted([int(x.split("_")[1]) for x in pdbs.keys()])	
	clusters_names = ["cluster_"+str(x) for x in cluster_rows]	

	with open('clusters.csv', 'wb') as csvfile:
		writer = csv.writer(csvfile, delimiter=';',quotechar='|', quoting=csv.QUOTE_MINIMAL)
		for cluster, c_name in zip(clusters_names, cluster_rows):
			cl_pdbs = ",".join(pdbs[cluster])
			cl_cogs = ",".join(cogs[cluster])
			cl_kogs = ",".join(kogs[cluster])
			cl_pfams = ",".join(pfams[cluster])
			writer.writerow([c_name, cl_pdbs, cl_pfams, cl_cogs, cl_kogs])

def get_cog_mapping(cog_file, kog_file, gis_dict):
	"""
	returns cog mappings for each cluster.
	cluster:[cogs]
	"""				
	db_file = open(CDD_DB, "r").readlines()
	cog = open(cog_mapping, "r").readlines()
	kog = open(kog_mapping, "r").readlines()
	cdd_db = {}

	for line in db_file:
		splitted = line.split()
		cdd_db[splitted[0]] = splitted[1]
	
	cog_dict = {}
	kog_dict = {}

	for ln in cog:
		cog_dict[ln.split()[0]] = ln.split()[1]
	
	for ln2 in kog:
		kog_dict[ln2.split()[0]] = ln2.split()[1]

	cluster2cog = {}
	cluster2kog = {}

	for cluster in gis_dict:
		cluster2cog[cluster] = []
		for gi in gi_dict[cluster]:
			if gi:
				if gi in cog_dict:
					cluster2cog[cluster].append(cdd_db[cog_dict[gi]])
		
		nonred = list(set(cluster2cog[cluster]))		
		cluster2cog[cluster] = nonred		
			
	for cluster in gis_dict:
		cluster2kog[cluster] = []
		for gi in gi_dict[cluster]:
			if gi:
				if gi in kog_dict:
					cluster2kog[cluster].append(cdd_db[kog_dict[gi]])
		
		nonred = list(set(cluster2kog[cluster]))	
		cluster2kog[cluster] = nonred

	return cluster2cog, cluster2kog		

def cluster2map(gis_dict):
	"""
	Make clusters mapping cytoscape file.
	"""
	outfile = open("clusters_mapping.txt", "w")
	
	for cluster in gis_dict:
		for gi in gis_dict[cluster]:
			outfile.write(gi+"\t"+cluster+"\n")

	outfile.close()		

def domains2graph_single(dom_dict):
	"""
	Prepare Cytoscape graph input for domains architecture.
	"""
	dir_path = os.path.join(final_cluster_dir, "graph_domains")

	if not os.path.exists(dir_path):
    		os.makedirs(dir_path)
	
	for cluster in dom_dict:
		cluster_combos = {}
		domains_file = open(os.path.join(final_cluster_dir,"graph_domains",cluster+".txt"), "w")
		for domain in dom_dict[cluster]:
			dom = list(set([x[3] for x in domain]))
			combos = [",".join(map(str,comb)) for comb in sorted(itertools.combinations(dom, 2))]
			for comb in combos:
				if comb not in cluster_combos:
					cluster_combos[comb]=0
				cluster_combos[comb]+=1
			
		for cl_combo in cluster_combos:
			domains_file.write("\t".join(cl_combo.split(","))+ "\t"+str(cluster_combos[cl_combo])+"\n")
		
		domains_file.close()

def domains2graph_all(dom_dict):
	domains_file = open(os.path.join(final_cluster_dir,"domains_graph.txt"), "w")
	
	cluster_combos = {}
	
	for cluster in dom_dict:
		for domain in dom_dict[cluster]:
			dom = list(set([x[3] for x in domain]))
			combos = [",".join(map(str,comb)) for comb in sorted(itertools.combinations(dom, 2))]
			for comb in combos:
				if comb not in cluster_combos:
					cluster_combos[comb]=0
				cluster_combos[comb]+=1
	for cl_combo in cluster_combos:
		domains_file.write("\t".join(cl_combo.split(","))+ "\t"+str(cluster_combos[cl_combo])+"\n")
	
	domains_file.close() 
		
def domains2GOmap(dom_dict, level):
	"""
	Map PFAM entries to EasyGo terms.
	"""
	
	gomap = open("easyGO_process_mapping.txt", "w")
	gos = open(EASY_GO, "r").read().splitlines()
	
	PFAMS = []
	for cluster in dom_dict:
		for dom in dom_dict[cluster]:
			PFAMS+=[x[3] for x in dom if x[3] != "His-Me"]
	
	nonred = list(set(PFAMS))
	
	
	GO_dict = {}
	
	for go in gos:
		splitted = go.split("\t")
		if level in splitted and splitted[0] in nonred:
			gomap.write(splitted[0]+"\t"+splitted[2]+"\n")
	
	gomap.close()

def domains2colors():
	"""
	Map PFAM domains to colors
	"""
	domains_list = open("pfam_domains.txt", "r").read().splitlines()
	colors =  random_colors.generate_N_pastels(len(domains_list))
	
	domains_col = {}
	
	for col, domain in zip(colors,domains_list):
		domains_col[domain] = col
	
	return domains_col

def functions2colors():
	"""
	Map functions to colors
	"""
	functions = open(FUNCTIONS, "r").read().splitlines()

	pfam2fun = {}
	fun2color = {}
	pfam2color = {}

	for fun in functions:
			splitted = fun.split(",")
			pfam2fun[splitted[0]] = splitted[1]

	values = list(set(pfam2fun.values()))
	colors =  random_colors.generate_N_pastels(len(values))
	
	for col, fun in zip(colors, values):
		fun2color[fun] = col


	for pfam in pfam2fun:
		pfam2color[pfam] = fun2color[pfam2fun[pfam]]

	get_colors_legend(fun2color)
	return pfam2color

def get_colors_legend(fun2color):
		"""
		Generate legend from colors and assigned functions.
		"""
		y = 0.5

		z = 2 #fixed width
		c = 0.6 #fixed height
		fig6 = plt.figure()
		ax6 = fig6.add_subplot(111, aspect='equal')

		for fun in fun2color:
			ax6.add_patch(patches.Rectangle((3,y),z,c,facecolor=fun2color[fun],linewidth=0))
		
			ax6.text(0.5*(3+10+z), 0.5*(y+y+c), fun,
					horizontalalignment='center', verticalalignment='center',
	        		fontsize=3, color='black')

			y+=0.8
		plt.ylim([0,35])
		plt.xlim([0,30])
		plt.axis('off')
		plt.savefig("function_legend.pdf", format='pdf')

		#fig6.savefig('legend.png', dpi=700, bbox_inches='tight')	

def get_protein_names(gi_dict, names_list):
	names = open(names_list, "r").read().splitlines()
	cl2protein = {}
	
	names_dict = {}
	for name in names:
		names_dict[name.split()[0]] = " ".join(name.split()[1:])
	
	for cluster in gi_dict:
		cl2protein[cluster] = []
		for gi in gi_dict[cluster]:
			if gi:
				try:
					cl2protein[cluster].append(names_dict[gi])
				except: pass
				
	for cluster in cl2protein:
		fn = open(cluster+"_protnames.txt", "w")
		for gi in cl2protein[cluster]:
			fn.write(gi+"\n")
		fn.close()		

def clusters2ws_proteins(gi_dict):
		
	all_seqs = {}	
	for rec in SeqIO.parse(open(whole_seqs), 'fasta'):
		all_seqs[rec.id] = str(rec.seq)
	
	for cluster in gi_dict:
		fh = open(cluster+"_ws.fas","w")
		for gi in gi_dict[cluster]:
			if gi in all_seqs:
				fh.write(">"+gi+"\n")
				fh.write(all_seqs[gi]+"\n")
		fh.close()

def clusters2taxonomy(gis_dict):
	taxfile = open("taxonomy_round_2.txt", "r").read().splitlines()
	tax_dict = {}
	
	for tax in taxfile:
		gi = tax.split()[0]
		taxa = tax.split()[1]
		tax_dict[gi] = taxa
	
	for cluster in gis_dict:
		fh = open(cluster+"_taxa.txt","w")
		list_tax = []
		for gi in gis_dict[cluster]:
			try:
				list_tax.append(tax_dict[gi])
			except: pass
		counts = Counter(list_tax)
		for count in counts:
			fh.write(count+"\t"+str(counts[count])+"\n")
		fh.close()

def cluster2organism2pfam(gis_dict):
	taxfile = open("gis2family_tax.txt","r").read().splitlines()
	pfamfile = open(pfam_mapping, "r").read().splitlines()
	
	taxdict = {}
	pfamdict = {}
	
	for tax in taxfile:
		splitted = tax.split()
		taxdict[splitted[0]] = splitted[1:]
	
	for pfam in pfamfile:
		splitted = pfam.split()
		pfamdict[splitted[0]] = splitted[1]
		
	for cluster in gis_dict:
		fh = open(cluster+"_tax_pfam.txt","w")
		for gi in gi_dict[cluster]:
			if gi in pfamdict and gi in taxdict:
				fh.write(gi+"\t"+pfamdict[gi]+"\t"+"\t".join(taxdict[gi])+"\n")
				
			elif gi in pfamdict and gi not in taxdict:
				fh.write(gi+"\t"+pfamdict[gi]+"\n")
			elif gi not in pfamdict and gi in taxdict:
				fh.write(gi+"\t"+"\t".join(taxdict[gi])+"\n")
			else: pass
		fh.close()

def do():	
	pass

if __name__ == "__main__":
	gi_dict = get_gis_dict(gis_dir)
	doms = domain_architectures(ws_pfam_mapping, gi_dict)
	#unfiltered_doms = domain_architectures(ws_pfam_mapping, gi_dict)
	#stats = domain2stats(unfiltered_doms)
	#get_protein_names(gi_dict, PROTEIN_NAMES)
	
	#clusters2ws_proteins(gi_dict)
	#clusters2taxonomy(gi_dict)
	#cluster2organism2pfam(gi_dict)
	#do rysowania najczestszych na posterze
	'''all = get_all_archs(unfiltered_doms)
	filtered = filter_archs(all)
	
	fl = sorted(filtered["all"], key=lambda x: x[-1])
	to_draw = {}
	to_draw["all"] = fl[-10:]
	draw_architectures(to_draw)'''
	
	#domains2file(unfiltered_doms)
	#doms2 = remove_repeats(doms)
	draw_architectures(doms)
	#get_taxonomy_count(taxonomy_file,gi_dict)
	#domains2graph_all(doms)
	#domains2GOmap(doms, "process")
	#all_gis = get_all_sequences(cd_hits, gi_dict)
	#pfams = get_pfam_mapping(gi_dict, pfam_mapping)	
	#pdbs = get_pdb90(all_gis, pdb_mappings)
	#cog, kog = get_cog_mapping(gi_dict, cog_mapping, gi_dict)
	#make_table(pdbs, pfams, cog, kog)

	#cluster2map(gi_dict)
	
	#filter_archs(doms)
	#map2limits(limits)
	#map2colors(gi_dict)	
	#gather_sequences(all_seqs_file, gi_dict)
	#run_mafft()

