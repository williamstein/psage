from sage.rings.all import NumberField, QQ

R = QQ['x']
x = R.gen()
F = NumberField(x*x - x- 1, 'a')
a = F.gen()
