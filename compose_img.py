# Script for composing the taxonomy image
# using PIL library. Custom lines spacing and 
# background size.

from PIL import Image
import glob
import ghostscript

def pdf2jpeg(pdf_input_path, jpeg_output_path):
    args = ["pdf2jpeg", # actual value doesn't matter
            "-dNOPAUSE",
            "-sDEVICE=jpeg",
            "-r144",
            "-sOutputFile=" + jpeg_output_path,
            pdf_input_path]
    ghostscript.Ghostscript(*args)

def glue_together(space=9):
	"""
	Compose taxonomy image. Indicate the spacing between the consecutive images.
	"""

	background = Image.new("RGB",(140, 842), (255, 255, 255, 255)) #A4
	o_w = 10 #offset values
	o_h = 20

	for i in range(1,52):
		img_name = "cluster_"+str(i)+"_taxa.png"
		img = Image.open(img_name, 'r') 
		img_w, img_h = img.size
		bg_w, bg_h = background.size
		offset = (o_w, o_h)
		o_h += space
		background.paste(img, offset)

	background.save('taxonomy_img.png', 'png')

if __name__ == "__main__":

	pdfs = glob.glob("*.png")

	'''for pdf in pdfs:
		out_name = pdf.split(".")[0]+".png"
		pdf2jpeg(pdf,out_name)'''

	glue_together()