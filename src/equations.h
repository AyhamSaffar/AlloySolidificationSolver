#ifndef EQUATIONS_H
#define EQUATIONS_H


namespace equation
{
    inline double getPt(double V, double R, double a);
    inline double getPc(double V, double R, double D);
    inline double getIvPt(double Pt, double E1Pt);
    inline double getIvPc(double Pc, double E1Pc);
    inline double getF1(double L, double Cp, double IvPt, double m, double C0, double k0, double IvPc, double r,
                        double R, double dT);
    inline double getQ(double Pt, double L, double Cp, double Pc, double m, double C0, double k0,double IvPc);
    inline double getF2(double r, double o, double Q, double R);
}

#endif
