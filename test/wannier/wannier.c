/* Analytical solution from:
   "A contribution to the hydrodynamics of lubrication"
   G.H. Wannier - Quarterly of Appl. Math., vol VIII, No.1, April 1950 */
void psiuv (double x, double y, 
	    double r1, double r2, double e,
	    double v1, double v2,
	    double * ux, double * uy, double * psi) {
  double d1 = (r2*r2 - r1*r1)/(2.*e) - e/2.;
  double d2 = d1 + e;
  double s = sqrt((r2 - r1 - e)*(r2 - r1 + e)*(r2 + r1 + e)*(r2 + r1 - e))/(2.*e);
  double l1 = log((d1 + s)/(d1 - s));
  double l2 = log((d2 + s)/(d2 - s));
  double den = (r2*r2 + r1*r1)*(l1 - l2) - 4.*s*e;
  double curlb = 2.*(d2*d2 - d1*d1)*(r1*v1 + r2*v2)/((r2*r2 + r1*r1)*den) +
    r1*r1*r2*r2*(v1/r1 - v2/r2)/(s*(r1*r1 + r2*r2)*(d2 - d1));
  double A = -0.5*(d1*d2-s*s)*curlb;
  double B = (d1 + s)*(d2 + s)*curlb;
  double C = (d1 - s)*(d2 - s)*curlb;
  double D = (d1*l2 - d2*l1)*(r1*v1 + r2*v2)/den
    - 2.*s*((r2*r2 - r1*r1)/(r2*r2 + r1*r1))*(r1*v1 + r2*v2)/den
    - r1*r1*r2*r2*(v1/r1 - v2/r2)/((r1*r1 + r2*r2)*e);
  double E = 0.5*(l1 - l2)*(r1*v1 + r2*v2)/den;
  double F = e*(r1*v1 + r2*v2)/den;

  y += d2;
  double spy = s + y;
  double smy = s - y;
  double zp = x*x + spy*spy;
  double zm = x*x + smy*smy;
  double l = log(zp/zm);
  double zr = 2.*(spy/zp + smy/zm);

  *psi = A*l + B*y*spy/zp + C*y*smy/zm + D*y + E*(x*x + y*y + s*s) + F*y*l;
  *ux = -A*zr - B*((s + 2.*y)*zp - 2.*spy*spy*y)/(zp*zp) -
    C*((s - 2.*y)*zm + 2.*smy*smy*y)/(zm*zm) - D -
    E*2.*y - F*(l + y*zr);
  *uy = -A*8.*s*x*y/(zp*zm) - B*2.*x*y*spy/(zp*zp) -
    C*2.*x*y*smy/(zm*zm) + E*2.*x -
    F*8.*s*x*y*y/(zp*zm);
}
