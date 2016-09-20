extern double nan_density(double t, double *pars, double *qp);
extern void nan_hessian(double t, double *pars, double *qp, double *hess);

extern double rotating_value(double t, double *Omega, double *qp);
extern void rotating_gradient(double t, double *Omega, double *qp, double *grad);

extern double henon_heiles_value(double t, double *pars, double *qp);
extern void henon_heiles_gradient(double t, double *pars, double *qp, double *grad);

extern double kepler_value(double t, double *pars, double *qp);
extern void kepler_gradient(double t, double *pars, double *qp, double *grad);
extern void kepler_hessian(double t, double *pars, double *qp, double *hess);

extern double isochrone_value(double t, double *pars, double *qp);
extern void isochrone_gradient(double t, double *pars, double *qp, double *grad);
extern double isochrone_density(double t, double *pars, double *qp);
extern void isochrone_hessian(double t, double *pars, double *qp, double *hess);

extern double hernquist_value(double t, double *pars, double *qp);
extern void hernquist_gradient(double t, double *pars, double *qp, double *grad);
extern double hernquist_density(double t, double *pars, double *qp);
extern void hernquist_hessian(double t, double *pars, double *qp, double *hess);

extern double plummer_value(double t, double *pars, double *qp);
extern void plummer_gradient(double t, double *pars, double *qp, double *grad);
extern double plummer_density(double t, double *pars, double *qp);
extern void plummer_hessian(double t, double *pars, double *qp, double *hess);

extern double jaffe_value(double t, double *pars, double *qp);
extern void jaffe_gradient(double t, double *pars, double *qp, double *grad);
extern double jaffe_density(double t, double *pars, double *qp);

extern double stone_value(double t, double *pars, double *qp);
extern void stone_gradient(double t, double *pars, double *qp, double *grad);
extern void stone_density(double t, double *pars, double *qp);

extern double sphericalnfw_value(double t, double *pars, double *qp);
extern void sphericalnfw_gradient(double t, double *pars, double *qp, double *grad);
extern double sphericalnfw_density(double t, double *pars, double *qp);
extern void sphericalnfw_hessian(double t, double *pars, double *qp, double *hess);

extern double flattenednfw_value(double t, double *pars, double *qp);
extern void flattenednfw_gradient(double t, double *pars, double *qp, double *grad);
extern double flattenednfw_density(double t, double *pars, double *qp);

extern double satoh_value(double t, double *pars, double *qp);
extern void satoh_gradient(double t, double *pars, double *qp, double *grad);
extern double satoh_density(double t, double *pars, double *qp);

extern double miyamotonagai_value(double t, double *pars, double *qp);
extern void miyamotonagai_gradient(double t, double *pars, double *qp, double *grad);
extern void miyamotonagai_hessian(double t, double *pars, double *qp, double *hess);
extern double miyamotonagai_density(double t, double *pars, double *qp);

extern double leesuto_value(double t, double *pars, double *qp);
extern void leesuto_gradient(double t, double *pars, double *qp, double *grad);
extern double leesuto_density(double t, double *pars, double *qp);

extern double logarithmic_value(double t, double *pars, double *qp);
extern void logarithmic_gradient(double t, double *pars, double *qp, double *grad);

extern double rotating_logarithmic_value(double t, double *pars, double *qp);
extern void rotating_logarithmic_gradient(double t, double *pars, double *qp, double *grad);

extern double lm10_value(double t, double *pars, double *qp);
extern void lm10_gradient(double t, double *pars, double *qp, double *grad);
