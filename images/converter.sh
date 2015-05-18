color=224,255,246

for image in theta.png prior.png posterior.png likelihood.png model.png bayes_rule.png dG.png dH.png H0.png Xs.png Mc.png sigma.png norm_n.png variance.png parameters.png
do
  convert ${image} -fuzz 100% -fill "rgb(${color})" -opaque black colored_${image}
done
