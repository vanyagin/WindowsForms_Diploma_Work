using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Timers;
using System.Threading;
using System.Windows.Forms.DataVisualization.Charting;

namespace WindowsForms_Diploma_Work
{




    public partial class Form1 : Form
    {




        class GausMethod
        {
            public uint RowCount;
            public uint ColumCount;
            public Complex[][] Matrix { get; set; }
            public Complex[] RightPart { get; set; }
            public Complex[] Answer { get; set; }

            public GausMethod(uint Row, uint Colum)
            {
                RightPart = new Complex[Row];
                Answer = new Complex[Row];
                Matrix = new Complex[Row][];
                for (int i = 0; i < Row; i++)
                    Matrix[i] = new Complex[Colum];
                RowCount = Row;
                ColumCount = Colum;

                //обнулим массив
                for (int i = 0; i < Row; i++)
                {
                    Answer[i] = 0;
                    RightPart[i] = 0;
                    for (int j = 0; j < Colum; j++)
                        Matrix[i][j] = 0;
                }
            }

            private void SortRows(int SortIndex)
            {

                Complex MaxElement = Matrix[SortIndex][SortIndex];
                int MaxElementIndex = SortIndex;
                for (int i = SortIndex + 1; i < RowCount; i++)
                {
                    if (Matrix[i][SortIndex].Magnitude > MaxElement.Magnitude)
                    {
                        MaxElement = Matrix[i][SortIndex];
                        MaxElementIndex = i;
                    }
                }

                //теперь найден максимальный элемент ставим его на верхнее место
                if (MaxElementIndex > SortIndex)//если это не первый элемент
                {
                    Complex Temp;

                    Temp = RightPart[MaxElementIndex];
                    RightPart[MaxElementIndex] = RightPart[SortIndex];
                    RightPart[SortIndex] = Temp;

                    for (int i = 0; i < ColumCount; i++)
                    {
                        Temp = Matrix[MaxElementIndex][i];
                        Matrix[MaxElementIndex][i] = Matrix[SortIndex][i];
                        Matrix[SortIndex][i] = Temp;
                    }
                }
            }

            public int SolveMatrix()
            {
                if (RowCount != ColumCount)
                    return 1; //нет решения

                for (int i = 0; i < RowCount - 1; i++)
                {
                    SortRows(i);
                    for (int j = i + 1; j < RowCount; j++)
                    {
                        if (Matrix[i][i] != 0) //если главный элемент не 0, то производим вычисления
                        {
                            Complex MultElement = Matrix[j][i] / Matrix[i][i];
                            for (int k = i; k < ColumCount; k++)
                                Matrix[j][k] -= Matrix[i][k] * MultElement;
                            RightPart[j] -= RightPart[i] * MultElement;
                        }
                        //для нулевого главного элемента просто пропускаем данный шаг
                    }
                }

                //ищем решение
                for (int i = (int)(RowCount - 1); i >= 0; i--)
                {
                    Answer[i] = RightPart[i];

                    for (int j = (int)(RowCount - 1); j > i; j--)
                        Answer[i] -= Matrix[i][j] * Answer[j];

                    if (Matrix[i][i] == 0)
                        if (RightPart[i] == 0)
                            return 2; //множество решений
                        else
                            return 1; //нет решения

                    Answer[i] /= Matrix[i][i];

                }
                return 0;
            }



            public override String ToString()
            {
                String S = "";
                for (int i = 0; i < RowCount; i++)
                {
                    S += "\r\n";
                    for (int j = 0; j < ColumCount; j++)
                    {
                        S += Matrix[i][j].ToString("F04") + "\t";
                    }

                    S += "\t" + Answer[i].ToString("F08");
                    S += "\t" + RightPart[i].ToString("F04");
                }
                return S;
            }
        }






        static double L_st = 0.02;
        static double l0_st = 0.01;
        static double theta_st = Math.PI / 4;
        static double d0_st = 0.5;
        static int k_st = 2;


        static double L;
        static double l0; // полудлина трещины (задается для частного случая)
        static double theta; //угол отклонения
        static double d0; //начальное смещение по оси x3
        const double v = 0.64; // c66 / c44
        const double h = 1; //толщина слоя

        static int k = 2;

        const int N = 20;
        const int M = N;
        const double ht = 0.1;
        const double p0 = 1;


        public static double tk(int k)
        {
            return -1 + (k - 1) * ht;
        }

        public static double tauk(int k)
        {
            return -1 + (k - 1) * ht + ht / 2;
        }

        public static double lambda0(int n)
        {
            return Math.PI * (n - 0.5) / h;
        }

        public static Complex lambda(int n)
        {
            return Math.PI * (n - 0.5) / h / Complex.ImaginaryOne;
        }

        public static Complex alpha(int n)
        {
            return Complex.Sqrt(k * k - lambda0(n) * lambda0(n)) / Math.Sqrt(v);
        }

        public static double q1(double t)
        {
            return l0 * Math.Cos(theta) * t;
        }

        public static double q3(double t)
        {
            return d0 + l0 * Math.Sin(theta) * t;
        }

        public static double q1_der()
        {
            return l0 * Math.Cos(theta);
        }

        public static double q3_der()
        {
            return l0 * Math.Sin(theta);
        }

        public static double c1()
        {
            return Math.Sqrt(v) * Math.Pow((v * q3_der() * q3_der() - q1_der() * q1_der()), 2) / (l0 * Math.Pow((v * q3_der() * q3_der() + q1_der() * q1_der()), 2));
        }

        public static double c2()
        {
            return 2 * v * Math.Sqrt(v) * (q3_der() * q1_der() + q1_der() * q3_der()) * q1_der() * q3_der() / (l0 * Math.Pow((v * q3_der() * q3_der() + q1_der() * q1_der()), 2));
        }

        public static double n1()
        {
            return q3_der() / l0;
        }

        public static double n3()
        {
            return q1_der() / l0;
        }

        public static Complex Fk(double tau)
        {
            int n = N;
            Complex sum = 0;
            for (int i = 1; i <= n; i++)
            {
                sum += Math.Pow(-1, i + 1) * (-Math.Sign(L + q1(tau)) * n1() * Math.Sin(lambda0(i) * q3(tau)) + (Complex.ImaginaryOne * lambda0(i) / v / alpha(i)) * n3() * Math.Cos(lambda0(i) * q3(tau))) * Complex.Pow(Math.E, Complex.ImaginaryOne * alpha(i) * Math.Abs(L + q1(tau)));
            }
            return sum;
        }

        public static Complex Aki(int k, int i)
        {
            return (c1() + c2()) * (1 / (tk(i) - tauk(k)) - 1 / (tk(i + 1) - tauk(k)));
            
        }

        public static Complex f1(int n)
        {
            return Complex.ImaginaryOne * v * q3_der() + lambda(n) * q1_der() / alpha(n);
        }

        public static Complex f2(int n)
        {
            return -Complex.ImaginaryOne * v * q3_der() + lambda(n) * q1_der() / alpha(n);
        }

        public static Complex[] Xsi = new Complex[10];

        static GausMethod Solution = new GausMethod(10, 10);

        static void SLAU_solve(double _L, double _l0, double _theta, double _d0)
        {
            L = _L;
            l0 = _l0;
            theta = _theta;
            d0 = _d0;

            //заполняем правую часть
            for (int i = 0; i < 10; i++)
            {
                Solution.RightPart[i] = Fk(tauk(i + 1));
            }

            //заполняем матрицу
            for (int i = 0; i < 10; i++)
                for (int j = 0; j < 10; j++)
                    Solution.Matrix[i][j] = Aki(i + 1, j + 1);

            //решаем матрицу
            Solution.SolveMatrix();

            //сохраняем ответ
            for (int i = 0; i < 10; i++)
            {
                Xsi[i] = Solution.Answer[i];

            }
        }

        static void SLAU_solve(double _L, double _l0, double _theta, double _d0, int _k)
        {
            k = _k;
            L = _L;
            l0 = _l0;
            theta = _theta;
            d0 = _d0;

            //заполняем правую часть
            for (int i = 0; i < 10; i++)
            {
                Solution.RightPart[i] = Fk(tauk(i + 1));
            }

            //заполняем матрицу
            for (int i = 0; i < 10; i++)
                for (int j = 0; j < 10; j++)
                    Solution.Matrix[i][j] = Aki(i + 1, j + 1);

            //решаем матрицу
            Solution.SolveMatrix();

            //сохраняем ответ
            for (int i = 0; i < 10; i++)
            {
                Xsi[i] = Solution.Answer[i];

            }
        }


        public static Complex An(int n)
        {


            Complex sum = 0;
            for (int i = 0; i < Xsi.Length; i++)
            {
                sum += Math.Pow(-1, n + 1) * Xsi[i] * (f1(n) * Complex.Exp(Complex.ImaginaryOne * lambda(n) * q3(tk(i))) + f2(n) * Complex.Exp(-Complex.ImaginaryOne * lambda(n) * q3(tk(i)))) * Complex.Exp(-Complex.ImaginaryOne * alpha(n) * q1(tk(i)));
            }
            return ht / 2 / h / v * sum;
        }


        public static Complex u(double x1)
        {
            Complex sum = 0;
            for (int i = 1; i <= 10; i++)
            {
                sum += An(i)* Complex.Exp(Complex.ImaginaryOne * alpha(i) * x1);
            }
            return sum;
        }


        public static double F(double _L, double _l0, double _theta, double _d0)
        {

            Complex sum = 0;
            for (int i = 1; i <= 2; i++)
            {
                SLAU_solve(L_st, l0_st, theta_st, d0_st);
                Complex _u_st = u(i);
                Console.WriteLine(_u_st);
                SLAU_solve(_L, _l0, _theta, _d0);
                Complex _u = u(i);
                Console.WriteLine(_u);
                sum += (_u_st - _u) * (_u_st - _u);
            }
            
            return sum.Magnitude;

        }

        public static double F1(double _L, double _l0, double _theta, double _d0)
        {

            double sum = 0;
            for (int i = 1; i <= 2; i++)
            {
                SLAU_solve(L_st, l0_st, theta_st, d0_st);
                Complex _u_st = u(i);
                //Console.WriteLine(_u_st);
                SLAU_solve(_L, _l0, _theta, _d0);
                Complex _u = u(i);
                //Console.WriteLine(_u);
                sum += Math.Abs((_u_st.Real - _u.Real)/_u_st.Real) + Math.Abs((_u_st.Imaginary - _u.Imaginary)/_u_st.Imaginary);
            }

            return sum;

        }


        public Form1()
        {
            InitializeComponent();


            


            Console.WriteLine(F1(L_st, l0_st, theta_st, d0_st));
            Console.WriteLine(F1(L_st, l0_st, theta_st, 0.9));






            SLAU_solve(L_st, l0_st, theta_st, d0_st, 1);
            //создаем элемент Chart
            Chart myChart = new Chart();
            //кладем его на форму и растягиваем на все окно.
            myChart.Parent = this;
            myChart.Dock = DockStyle.Fill;
            //добавляем в Chart область для рисования графиков, их может быть
            //много, поэтому даем ей имя.
            myChart.ChartAreas.Add(new ChartArea("Math functions"));
            //Создаем и настраиваем набор точек для рисования графика, в том
            //не забыв указать имя области на которой хотим отобразить этот
            //набор точек.
            Series s1 = new Series("1");
            Series s2 = new Series("2");
            Series s3 = new Series("3");
            s1.ChartType = SeriesChartType.Line;
            s1.ChartArea = "Math functions";
            s1.Color = Color.Red;
            s2.ChartType = SeriesChartType.Line;
            s2.ChartArea = "Math functions";
            s2.Color = Color.Green;
            s3.ChartType = SeriesChartType.Line;
            s3.ChartArea = "Math functions";
            s3.Color = Color.Blue;
            for (double i = 0.3; i <= 10; i += 0.1)
            {
                s1.Points.AddXY(i, u(i).Real);
            }
            myChart.Series.Add(s1);


            SLAU_solve(L_st, l0_st, theta_st, d0_st, k_st);


            for (double i = 0.3; i <= 10; i += 0.1)
            {
                s2.Points.AddXY(i, u(i).Real);
            }
            myChart.Series.Add(s2);


            SLAU_solve(L_st, l0_st, theta_st, d0_st, 6);


            for (double i = 0.3; i <= 10; i += 0.1)
            {
                s3.Points.AddXY(i, u(i).Real);
            }
            myChart.Series.Add(s3);
        }

        private void button1_Click(object sender, EventArgs e)
        {
            //visible = false
        }


    }

    
}
