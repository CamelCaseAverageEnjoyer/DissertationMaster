from PIL import Image

for i in range(1757):  # 1757
    im = Image.open(f"img/gibbon_{(i+1):04}.png")
    im1 = im.crop((600, 390, 1370, 760))
    im1.save(f"img2/gibbon_{(i+1):04}.png")