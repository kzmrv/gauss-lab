using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Providers;
using MathNet.Numerics.LinearAlgebra.Solvers;

namespace lab02
{
    class Program
    {
        public static void WriteMatrix<T>(Matrix<T> matrix, string name = "") where T : struct,IEquatable<T>, IFormattable
        {
            Console.WriteLine("Output matrix {0}", name);
            for (int j = 0; j < matrix.ColumnCount; j++)
            {
                foreach (T el in matrix.Row(j))
                {
                    Console.Write(string.Format("{0,10:N5}", el));
                }
                Console.Write('\n');
            }
            Console.WriteLine("**********************");
        }

        public static void WriteVector<T>(Vector<T> vector, string name) where T : struct,IEquatable<T>, IFormattable
        {
            Console.WriteLine("Output vector : {0}", name);
            for (int i = 0; i < vector.Count; i++)
            {
                Console.Write(string.Format("{0,20:N5}", vector[i]));
            }
            Console.WriteLine();
            Console.WriteLine("**********************");

        }
        public static void Gauss(Matrix<double> A, Vector<double> B, out Vector<double> X)
        {
            int n = A.RowCount;
            var inputA = Matrix<double>.Build.Dense(n, n);
            var inputB = Vector<double>.Build.Dense(n);
            A.CopyTo(inputA);
            B.CopyTo(inputB);
            WriteMatrix(A, "A");
            WriteVector(B, "B");
            Console.Write("\n\n");
            for (int i = 0; i < n; i++)
            {
                Console.WriteLine(" *** Iteration {0}", i + 1);
                WriteMatrix(A, "A" + i.ToString());
                int maxIndex = A.Diagonal().SubVector(i, n - i).AbsoluteMaximumIndex() + i;
                var P = Matrix<double>.Build.DenseIdentity(n, n);
                P[maxIndex, maxIndex] = 0;
                P[i, i] = 0;
                P[maxIndex, i] = 1;
                P[i, maxIndex] = 1;
                WriteMatrix(P, "P" + i.ToString());
                A = P * A;
                WriteMatrix(A, "AT" + i.ToString());
                Console.WriteLine(maxIndex);
                var M = Matrix<double>.Build.DenseIdentity(n, n);
                M[i, i] = 1 / A[i, i];
                for (int j = i + 1; j < n; j++)
                {
                    M[j, i] = -A[j, i] / A[i, i];
                }
                WriteMatrix(M, "M" + i.ToString());
                A = M * A;
                B = P * B;
                B = M * B;
                Console.WriteLine("B" + i.ToString());
                Console.WriteLine(B);
                Console.Write("\n\n");
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
            Console.WriteLine("X");
            Console.WriteLine(X);
            var result = inputA * X - inputB;
            WriteVector(result, "A*X - B");
        }

        public static void Jacobi(Matrix<double> A, Vector<double> B, out Vector<double> X)
        {
            double eps = 0.0001;
            int n = A.RowCount;
            Matrix<double> DR = (Matrix<double>.Build.Diagonal(A.Diagonal().ToArray())).Inverse();
            WriteMatrix(DR, "Diagonal reverse");
            Matrix<double> U = A.StrictlyUpperTriangle();
            WriteMatrix(U, "Upper triangle");
            Matrix<double> L = A.StrictlyLowerTriangle();
            WriteMatrix(L, "Lower triangle");
            Matrix<double> SP = -DR * (L + U);
            Vector<double> prevX = Vector<double>.Build.Dense(n);
            X = Vector<double>.Build.Dense(new double[] { 8e6, 6e6, -3e6, 2e6 });
            do
            {
                X.CopyTo(prevX);
                X = SP * prevX -DR * B;
            }
            while ((prevX - X).L2Norm() > eps);



        }
        static void Main(string[] args)
        {
            int n = 4;
            var B = Vector<double>.Build.Dense(n, 24);
            var C = Matrix<double>.Build.Dense(n, n, (i, j) => 2.4 * (i + 1) * (j + 1));
            var A = Matrix<double>.Build.Dense(n, n, (i, j) => 159 / (10 * C[i, j] * C[i, j] *
                C[i, j] + C[i, j] * C[i, j] + 25));
            var inputA = Matrix<double>.Build.Dense(n, n);
            A.CopyTo(inputA);
            Vector<double> gaussRes;
            Vector<double> jacobiRes;
            Gauss(A, B, out gaussRes);
            WriteVector(gaussRes, "Gauss Result");
            WriteMatrix(A);
            Jacobi(A, B, out jacobiRes);
            WriteVector(jacobiRes, "Jacobi result");
            Console.Read();
            

        }
    }
}
