color=255,255,255
for image in latex_images/* 
do
convert ${image} -fuzz 100% -fill "rgb(${color})" -opaque black colored_${image##*/}
done
