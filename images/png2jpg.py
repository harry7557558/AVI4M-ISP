import os
from PIL import Image

files = [s for s in os.listdir('original/') if s[-4:]=='.png']

for file in files:
    img = Image.open('original/'+file).convert('RGB')
    if ' ' in file:
        file = file.split()[1]
    file = file.lower().replace('_', '-')
    print(file)
    file_png = 'png/'+file
    file_jpg = 'jpg/'+file.replace('.png', '.jpg')
    img.save(file_jpg)
