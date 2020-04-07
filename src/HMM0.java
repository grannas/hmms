import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

public class HMM0 {
    public static void main(String[] args) throws IOException {

        BufferedReader bReader = new BufferedReader(new InputStreamReader(System.in));

        // Kattis input
        double[][] A = inputReader(bReader);
        double[][] B = inputReader(bReader);
        double[][] pi = inputReader(bReader);

        // Kattis calculations
        double[][] A1 = matrixVectorMultiplication(pi, A);
        double[][] A2 = matrixVectorMultiplication(A1, B);

        // Kattis output
        System.out.print(A2.length + " ");
        System.out.print(A2[0].length + " ");
        for (int i = 0; i < A2[0].length; i++) {
            System.out.print(A2[0][i] + " ");
        }
    }

    private static double[][] matrixVectorMultiplication(double[][] vector, double[][] matrix) {

        int vectorCols = vector[0].length;
        int matrixCols = matrix[0].length;
        int matrixRows = matrix.length;

        double[][] newVector = new double[1][matrixCols];

        if (vectorCols != matrixRows) {
            System.err.println("The vector and the matrix dimensions do not match for multiplication");
        }  else {
            for (int k = 0; k < matrixCols; k++) {
                double sum = 0;

                for (int j = 0; j < matrixRows; j++) {
                    sum += vector[0][j] * matrix[j][k];
                }
                newVector[0][k] = sum;
            }
        }
        return newVector;
    }

    private static double[][] inputReader(BufferedReader bReader) throws IOException {
        System.err.print("Input matrix: ");
        String line = bReader.readLine();
        String[] inputArray = line.split(" ");

        int rows = Integer.parseInt(inputArray[0]);
        int cols = Integer.parseInt(inputArray[1]);

        double[][] matrix = new double[rows][cols];

        // Index in original array
        int index = 2;

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                matrix[i][j] = Double.parseDouble(inputArray[index]);
                index++;
            }
        }
        return matrix;
    }
}
