using System;
using System.Numerics;
using Accord.Math.Integration;

namespace SpectrumHHG
{
    public class Program
    {
        const double Ip = 21.55;
        const double Up = 57.38;
        const double C = -0.1043234466544098;
        const double Tir = 20;
        const double Txuv = 0.55;
        const double Wxuv = 30;

        double Fir0 => 2 * sqrt(Up);
        double Fxuv0 => 2 * Wxuv * sqrt(0.006273492); // Up_xuv = 0.006273492

        double sin(double t) => Math.Sin(t);
        double cos(double t) => Math.Cos(t);
        double sqrt(double t) => Math.Sqrt(t);
        double sqr(double t) => Math.Pow(t, 2);
        double pow(double x, double y) => Math.Pow(x, y);
        double exp(double t) => Math.Exp(t);
        double ln(double t) => Math.Log(t);
        double PI => Math.PI;

        Complex a_prop(double OMEGA, double t1, double t2)
        {
            var exp = Complex.Exp(new Complex(0, -S(t1, t2) + (OMEGA * t2)));
            return exp / sqrt(pow(t2 - t1, 3));
        }

        Complex a_ion(double t1, double t2)
        {
            var exp = Complex.Exp(new Complex(0, -Wxuv * t1));
            double a1 = K(t1, t1, t2) / (t2 - t1);
            double a2 = 2 * PI * sqrt(K(t1, t1, t2) * (F(t1) + a1));
            return -d(t1, t2) * exp / a2;
        }

        Complex a_rec(double t1, double t2)
        {
            return Complex.Conjugate(d(t1, t2));
        }

        double S(double t1, double t2)
        {
            double s1 = Ip * (t2 - t1);
            double s2 = TrapezoidalRule.Integrate(t => sqr(K(t, t1, t2)), t1, t2, (int)Math.Round(200 * (t2 - t1)));
            return s1 + (s2 / 2);
        }

        double K(double t, double t1, double t2)
        {
            return A(t) + Q(t1, t2);
        }

        double Q(double t1, double t2)
        {
            if (t1 == t2)
                return 0;

            return -(IntA(t2) - IntA(t1)) / (t2 - t1);
        }

        double Fxuv(double t, double t1)
        {
            t -= t1;
            return Fxuv0 * cos(Wxuv * t) * exp(-2 * ln(2) * sqr(t) / sqr(Txuv));
        }

        double F(double t)
        {
            if (t < 0 || t > Tir)
                return 0;

            return Fir0 * cos(t) * pow(sin(PI * t / Tir), 2);
        }

        double A(double t)
        {
            if (t < 0 || t > Tir)
                return 0;

            double a1 = (2 * PI / Tir) + 1;
            double a2 = (2 * PI / Tir) - 1;

            double y1 = -sin(t) / 2;
            double y2 = sin(a1 * t) / (4 * a1);
            double y3 = sin(a2 * t) / (4 * a2);
            return Fir0 * (y1 + y2 + y3 + C);
        }

        double IntA(double t)
        {
            if (t < 0 || t > Tir)
                return 0;

            double a1 = (2 * PI / Tir) + 1;
            double a2 = (2 * PI / Tir) - 1;

            double y1 = cos(t) / 2;
            double y2 = -cos(a1 * t) / (4 * a1 * a1);
            double y3 = -cos(a2 * t) / (4 * a2 * a2);
            return Fir0 * (y1 + y2 + y3 + (C * t));
        }

        Complex d(double t1, double t2)
        {
            return Complex.One;
        }

        public static void Main(string[] args)
        {
            Console.ReadLine();
        }
    }
}