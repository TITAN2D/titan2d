from PIL import Image
from PIL import ImageDraw 
from PIL import ImageFont
from os import chdir, path
import sys

def makeLabel(blankImage, fontDir, outImage, maxHeight,bg="#000000",fg="#ffffff",font="FreeSerif.ttf",FontSize=12):
    font_dir = fontDir
    img_name = outImage
    font_size = FontSize
    font_path = font_dir+font
    fnt = ImageFont.truetype(font_path, font_size)
    lineWidth = 20
    img = Image.open(blankImage)
    draw = ImageDraw.Draw(img)                     # setup to draw on the main image
    for x in range(11):
        text = "%.2f" % (maxHeight/10 *(10-x))
        draw.text((10, 33 +22*x), text, font=fnt, fill="black")      # add some text to the main
    del draw 
    img.save(img_name,"JPEG",quality=100)  
