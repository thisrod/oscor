r := 8/5;
u := 0.75*2pc;  h := 0.75*7pc;  l := -0.75*4pc;

color gris;  gris = 0.6 white;
transform sq, qs;
(0,0) transformed sq = (0,0) transformed qs = (0,0);
(0,1) transformed sq = (0,r);  (r,0) transformed sq = (1,0);
(0,r) transformed qs = (0,1);  (1,0) transformed qs = (r,0);

beginfig(1)

for t = 0 upto 6:
	fill fullcircle transformed sq rotated 30t transformed qs scaled u shifted (t*h,0) withcolor gris;
	fill fullcircle transformed sq rotated 30t scaled u shifted (t*h,l) withcolor gris;
endfor

label.lft(btex Bose field etex, (-u/2-1pc,0));
label.lft(btex sound wave etex, (-u/2-1pc,l));
x1 = 0;  y1 = l-r*u/2-0.5pc;  z2 = z1 + 6h*right;
pair ofst;  ofst = 3/4u*right;
drawarrow (z1-ofst)--(z2+ofst);
label.rt(btex $t$ etex, z2+ofst);
dotlabel.bot("0", z1);  dotlabel.bot(btex $2\pi\xi/\kappa$ etex, z2);

endfig
end