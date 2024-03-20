"""Я создан потому что хозяин почему-то решил что сохранять в png нормально"""
from PIL import Image

# Уменьшить число кадров
n = 8614
rate = 6
for i in range(int(n // rate)):
        im = Image.open(f"img/10/gibbon_{((i + 1) * rate):0{4}}.png") 
        im.save(f"img2/10/gibbon_{(i + 1):0{4}}.jpg")
