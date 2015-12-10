using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using MathNet.Numerics.LinearAlgebra;

namespace lab02
{
    class Program
    {
        public static void WriteLine(string s) {

        }
        public static void WriteMatrix<T>(Matrix<T> matrix, StreamWriter sw, string name = "") where T : struct,IEquatable<T>, IFormattable
        {
            sw.WriteLine("Output matrix {0}", name);
            for (int j = 0; j < matrix.ColumnCount; j++)
            {
                foreach (T el in matrix.Row(j))
                {
                    sw.Write(string.Format("{0,14:N8}", el));
                }
                sw.WriteLine();
            }
            sw.WriteLine("**********************");
        }

        public static void WriteVector<T>(Vector<T> vector, string name, StreamWriter sw) where T : struct,IEquatable<T>, IFormattable
        {

            sw.WriteLine("Output vector : {0}", name);
            for (int i = 0; i < vector.Count; i++)
            {
                sw.Write(string.Format("{0,20:N16}", vector[i]));
            }
            sw.WriteLine();
            sw.WriteLine("**********************");

        }
        public static void Gauss(Matrix<double> A, Vector<double> B, out Vector<double> X)
        {
            int n = A.RowCount;
            var inputA = Matrix<double>.Build.Dense(n, n);
            var inputB = Vector<double>.Build.Dense(n);
            A.CopyTo(inputA);
            B.CopyTo(inputB);
            for (int i = 0; i < n; i++)
            {
                int maxIndex = A.Diagonal().SubVector(i, n - i).AbsoluteMaximumIndex() + i;
                var P = Matrix<double>.Build.DenseIdentity(n, n);
                P[maxIndex, maxIndex] = 0;
                P[i, i] = 0;
                P[maxIndex, i] = 1;
                P[i, maxIndex] = 1;
                A = P * A;
                var M = Matrix<double>.Build.DenseIdentity(n, n);
                M[i, i] = 1 / A[i, i];
                for (int j = i + 1; j < n; j++)
                {
                    M[j, i] = -A[j, i] / A[i, i];
                }

                A = M * A;
                B = P * B;
                B = M * B;
                if (i == n - 2)
                {
                    double det = 1;
                    foreach (double el in A.Diagonal())
                    {
                        det *= el;
                    }
                    WriteLine("Determinant Gauss }"+det.ToString());
                }

            }

            X = Vector<double>.Build.Dense(n);
            for (int i = n - 1; i >= 0; i--)
            {
                double res = B[i];
                for (int j = n - 1; j > i; j--)
                {
                    res -= A[i, j] * X[j];
                }
                X[i] = res;
            }
        }
        public static void Jacobi(Matrix<double> A, Vector<double> B, out Vector<double> X)
        {
            double eps = 0.0001;
            int n = A.RowCount;
            Matrix<double> DR = (Matrix<double>.Build.Diagonal(A.Diagonal().ToArray())).Inverse();
            Matrix<double> U = A.StrictlyUpperTriangle();
            Matrix<double> L = A.StrictlyLowerTriangle();
            Matrix<double> SP = -DR * (L + U);
            Vector<double> prevX = Vector<double>.Build.Dense(n);
            X = Vector<double>.Build.Dense(new double[] { 8e6, 6e6, -3e6, 2e6, 2e6, 2e6, 2e6 });
            do
            {
                X.CopyTo(prevX);
                X = SP * prevX - DR * B;
            }
            while ((prevX - X).L2Norm() > eps);
            X = -X;
        }
        static void Main(string[] args)
        {
            StreamWriter sw = new StreamWriter(@"c:\Vasili\output.txt");
            int n = 7;
            var B = Vector<double>.Build.Dense(n, 24);
            var C = Matrix<double>.Build.Dense(n, n, (i, j) => 2.4 * (i + 1) * (j + 1));
            var A = Matrix<double>.Build.Dense(n, n, (i, j) => 159 / (10 * C[i, j] * C[i, j] *
                C[i, j] + C[i, j] * C[i, j] + 25));
            for (int i = 0; i < n; i++)
            {
                A[i, i] += 10;
            }
            WriteMatrix(A, sw, "A");
            WriteVector(B, "B", sw);
            var inputA = Matrix<double>.Build.Dense(n, n);
            A.CopyTo(inputA);
            Vector<double> gaussRes;
            Vector<double> jacobiRes;
            Gauss(A, B, out gaussRes);
            WriteVector(gaussRes, "Gauss Result", sw);
            Console.WriteLine("Determinant {0} ", A.Determinant());
            Jacobi(A, B, out jacobiRes);
            WriteVector(jacobiRes, "Jacobi result", sw);
            WriteVector(B - A * gaussRes, " diff Gauss", sw);
            WriteVector(B - A * jacobiRes, "diff Jacobi", sw);
            WriteMatrix(A.Inverse(), sw, "Inverse A");

            var Anew = Matrix<double>.Build.Dense(n, n);
            A.CopyTo(Anew);
            var AR = A.Inverse();
            double sumA = 0;
            double sumAR = 0;
            for (int i = 0; i < n; i++)
            {
                sumA += A.Row(i).Maximum();
                sumAR += AR.Row(i).Maximum();
            }
            sw.WriteLine("Cond : {0} ", sumA * sumAR);
            double eps = 0.000001;
            for (int i = 0; i < n; i++)
            {
                Anew[1, i] += eps;
            }
            var jacobiRes2 = Vector<double>.Build.Dense(n);
            var gaussRes2 = Vector<double>.Build.Dense(n);
            WriteMatrix(AR * A, sw, "A-1 * A");
            Gauss(Anew, B, out gaussRes2);
            Jacobi(Anew, B, out jacobiRes2);
            WriteVector(gaussRes2, " X* ", sw);
            double err = (gaussRes2 - gaussRes).L2Norm() / gaussRes.L2Norm();
            double err2 = (jacobiRes2 - jacobiRes).L2Norm() / jacobiRes.L2Norm();
            sw.WriteLine("Gauss error {0}", err);
            sw.WriteLine("Jacobi error {0}", err2);
            sw.Close();
            Console.WriteLine("Press any key to finish ...");
            Console.Read();
        }
    }
}
