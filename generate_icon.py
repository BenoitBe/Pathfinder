import os
import numpy as np
from scipy.ndimage import gaussian_filter
from PIL import Image, ImageDraw, ImageFont, ImageFilter

SS = 4
W, H = 200 * SS, 240 * SS
rng = np.random.default_rng(11)
SERIF = "/usr/share/fonts/truetype/liberation/LiberationSerif-Regular.ttf"

def P(fx, fy):
    return (fx * W, fy * H)

def fractal(h, w, sigmas=(2, 5, 12, 26), amps=(0.5, 0.3, 0.15, 0.08)):
    acc = np.zeros((h, w))
    for s, a in zip(sigmas, amps):
        acc += gaussian_filter(rng.standard_normal((h, w)), sigma=s * SS / 4) * a
    return (acc - acc.min()) / (np.ptp(acc) + 1e-9)

yy, xx = np.mgrid[0:H, 0:W].astype(float)
rad = np.sqrt((xx - W*0.5)**2 + (yy - H*0.42)**2) / (W*0.85)
img = np.array([5,7,15]) + (np.array([14,21,38]) - np.array([5,7,15])) * np.clip(1-rad,0,1)[...,None]
hy = H*0.60 + (xx - W*0.5)**2 / (W*4.5)
ground = yy >= hy
star = np.zeros((H, W)); ns = 600
sx, sy = rng.integers(0, W, ns), rng.integers(0, int(H*0.60), ns)
star[sy, sx] = rng.uniform(0.2, 1.0, ns)
star = gaussian_filter(star, sigma=0.6*SS)
img += star[...,None] * np.array([150,180,220]) * 0.6
gnz = fractal(H, W)
depth = np.clip((yy - hy) / (H*0.40), 0, 1)
gcol = (np.array([168,96,56]) * (1-depth[...,None]) + np.array([86,44,26]) * depth[...,None])
gcol = gcol * (0.8 + 0.4*gnz[...,None])
img = np.where(ground[...,None], gcol, img)
haze = np.exp(-((yy - hy)/(H*0.03))**2) * 0.5
img += (haze)[...,None] * np.array([180,120,110]) * 0.6
img += (fractal(H, W, (1.2,), (1.0,)) - 0.5)[...,None] * 6
img = np.clip(img, 0, 255).astype(np.uint8)

pil = Image.fromarray(img, "RGB").convert("RGBA")

shadow = Image.new("RGBA", (W, H), (0,0,0,0))
ImageDraw.Draw(shadow).ellipse([W*0.30, H*0.60, W*0.72, H*0.645], fill=(0,0,0,120))
shadow = shadow.filter(ImageFilter.GaussianBlur(3*SS))
pil = Image.alpha_composite(pil, shadow)
d = ImageDraw.Draw(pil, "RGBA")

GOLD, GOLD_HI, GOLD_SH = (200,152,72,255), (236,200,126,255), (120,84,36,255)
PANEL, PANEL_LN, PANEL_HI = (42,56,84,255), (22,30,48,255), (110,140,178,255)
WHEEL, WHEEL_D, HUB = (66,66,72,255), (28,28,34,255), (120,120,128,255)
STRUT = (170,170,178,255)

def wheel(cx, cy, r):
    d.ellipse([cx-r, cy-r*0.92, cx+r, cy+r*0.92], fill=WHEEL, outline=WHEEL_D, width=SS)
    for k in range(12):
        a = k/12*2*np.pi
        d.line([(cx+np.cos(a)*r*0.6, cy+np.sin(a)*r*0.85*0.6),
                (cx+np.cos(a)*r*0.95, cy+np.sin(a)*r*0.85*0.95)], fill=WHEEL_D, width=SS)
    d.ellipse([cx-r*0.28, cy-r*0.28, cx+r*0.28, cy+r*0.28], fill=HUB)
    d.arc([cx-r, cy-r*0.92, cx+r, cy+r*0.92], 200, 320, fill=(150,150,158,200), width=SS)

wy = 0.585; wxs = [0.375, 0.505, 0.635]; wr = W*0.055
for wx in [0.44, 0.57]:
    d.ellipse([P(wx-0.03, 0)[0], P(0,wy-0.05)[1], P(wx+0.03,0)[0], P(0,wy+0.028)[1]],
              fill=(50,50,56,255), outline=WHEEL_D, width=SS)
pivot = P(0.52, 0.50)
for wx in wxs:
    d.line([P(wx, wy), pivot], fill=STRUT, width=int(2.4*SS))
d.line([P(0.505, wy), P(0.635, wy)], fill=STRUT, width=int(2.4*SS))
for wx in wxs:
    wheel(*P(wx, wy), wr)

d.polygon([P(0.375,0.475), P(0.635,0.475), P(0.615,0.55), P(0.395,0.55)], fill=GOLD)
d.line([P(0.375,0.475), P(0.635,0.475)], fill=GOLD_HI, width=int(1.6*SS))
d.line([P(0.395,0.55), P(0.615,0.55)], fill=GOLD_SH, width=int(1.6*SS))
for fx in np.linspace(0.41, 0.60, 4):
    d.line([P(fx,0.48), P(fx-0.01,0.545)], fill=(230,196,126,90), width=SS)

d.polygon([P(0.345,0.455), P(0.675,0.432), P(0.665,0.462), P(0.335,0.485)], fill=PANEL)
d.polygon([P(0.335,0.485), P(0.665,0.462), P(0.665,0.472), P(0.335,0.495)], fill=(26,34,52,255))
for t in np.linspace(0.08, 0.92, 7):
    d.line([ (P(0.345,0.455)[0]*(1-t)+P(0.675,0.432)[0]*t, P(0.345,0.455)[1]*(1-t)+P(0.675,0.432)[1]*t),
             (P(0.335,0.485)[0]*(1-t)+P(0.665,0.462)[0]*t, P(0.335,0.485)[1]*(1-t)+P(0.665,0.462)[1]*t)],
           fill=PANEL_LN, width=max(1,SS//2))
d.line([P(0.345,0.470), P(0.675,0.447)], fill=PANEL_LN, width=max(1,SS//2))
d.line([P(0.40,0.462), P(0.62,0.442)], fill=PANEL_HI, width=SS)

d.line([P(0.41,0.475), P(0.40,0.37)], fill=STRUT, width=int(1.6*SS))
tip = P(0.40,0.365)
d.ellipse([tip[0]-2*SS, tip[1]-2*SS, tip[0]+2*SS, tip[1]+2*SS], fill=(220,220,228,255))

def ticks(y0, n=46, region=(0.10,0.90), hmax=16, up=True):
    x0, x1 = W*region[0], W*region[1]
    for i in range(n):
        x = x0 + (x1-x0)*i/(n-1)
        hh = rng.uniform(0.15,1.0)**1.5 * hmax * SS
        col = (127,212,255,255) if rng.random()>0.82 else (47,109,148,210)
        yb = y0-hh if up else y0+hh
        d.line([(x,y0),(x,yb)], fill=col, width=max(1,SS//2))
ticks(H*0.135, up=True); ticks(H*0.135+3*SS, up=False, hmax=6)
ticks(H*0.865, up=False); ticks(H*0.865-3*SS, up=True, hmax=6)

def spaced(cx, cy, text, font, sp, fill):
    ws=[d.textlength(c,font=font) for c in text]; tot=sum(ws)+sp*(len(text)-1); x=cx-tot/2
    for c,wd in zip(text,ws): d.text((x,cy),c,font=font,fill=fill); x+=wd+sp
title="PATHFINDER"; fs=int(20*SS); font=ImageFont.truetype(SERIF,fs); sp=int(3.5*SS)
while sum(d.textlength(c,font=font) for c in title)+sp*(len(title)-1) > W*0.80 and fs>10*SS:
    fs-=SS; font=ImageFont.truetype(SERIF,fs)
spaced(W/2, H*0.885, title, font, sp, (200,160,112,255))
spaced(W/2, H*0.855, "· MARS · PATHFINDER · SOJOURNER ·",
       ImageFont.truetype(SERIF,int(5.2*SS)), int(1.6*SS), (120,140,160,200))

mask=Image.new("L",(W,H),0)
ImageDraw.Draw(mask).rounded_rectangle([2*SS,2*SS,W-2*SS,H-2*SS], radius=16*SS, fill=255)
frame=Image.new("RGBA",(W,H),(0,0,0,0))
ImageDraw.Draw(frame).rounded_rectangle([2*SS,2*SS,W-2*SS,H-2*SS], radius=16*SS,
                                        outline=(60,96,130,180), width=max(1,SS))
pil=Image.composite(pil, Image.new("RGBA",(W,H),(0,0,0,0)), mask)
pil=Image.alpha_composite(pil, frame)
os.makedirs("assets", exist_ok=True)
pil.save("assets/pathfinder_icon_800.png")
for size in (512,256,200,128):
    pil.resize((int(size*W/H), size), Image.LANCZOS).save(f"assets/pathfinder_icon_{size}.png")
pil.resize((205,246), Image.LANCZOS).save("assets/pathfinder_icon.ico",
                                          sizes=[(64,77),(128,154),(205,246)])
print("Rover icon generated -> assets/")
