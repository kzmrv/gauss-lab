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
        public static void WriteMatrix(Matrix<double> matrix, string name="")
        {
            Console.WriteLine("Output matrix {0}", name);
            for (int j = 0; j < matrix.ColumnCount; j++)
            {
                foreach (double el in matrix.Row(j))
                {
                    Console.Write(string.Format("{0,10:N5}", el));
                }
                Console.Write('\n');
            }
            Console.WriteLine("**********************");
        }

        public static void WriteVector(Vector<double> vector, string name)
        {
            Console.WriteLine("Output vector : {0}", name);
            for (int i = 0; i < vector.Count; i++)
            {
                Console.Write(string.Format("{0,10:N5}", vector[i]));
            }
            Console.WriteLine();
            Console.WriteLine("**********************");

        }

        static void Main(string[] args)
        {
            int n = 4;
            var B = Vector<double>.Build.Dense(n, 24);
            var C = Matrix<double>.Build.Dense(n, n, (i, j) => 2.4 * (i+1) * (j+1));
            var A = Matrix<double>.Build.Dense(n, n, (i, j) => 159 / (10 * C[i, j] * C[i, j] *
                C[i, j] + C[i, j] * C[i, j] + 25));
            var inputA = Matrix<double>.Build.Dense(n,n);
            A.CopyTo(inputA);
            WriteVector(B, "B");
            WriteMatrix(A, "A");
            WriteMatrix(C, "C");
            /*
            A.SetColumn(0, new double[] { 10, 3, -2 });
            A.SetColumn(1, new double[] { 0, -1, 4 });
            A.SetColumn(2, new double[] { 3, 0, 1 });
            B = Vector<double>.Build.Dense(new double[] { 7, 2, 1 });
            */
            for (int i = 0; i < n; i++)
            {
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
            }
            var X = Vector<double>.Build.Dense(n);
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
            var result = inputA * X;
            WriteVector(result, "Res");
            Console.Read();
        }
    }
}
