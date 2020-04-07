import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;

public class HMM1 {
    public static void main(String[] args) throws IOException {

        BufferedReader bReader = new BufferedReader(new InputStreamReader(System.in));

        // Kattis input
        double[][] A = matrixInputReader(bReader);
        double[][] B = matrixInputReader(bReader);
        double[][] pi = matrixInputReader(bReader);

        String[] oTemp= bReader.readLine().split(" ");
        int observations = Integer.parseInt(oTemp[0]);

        int[] O = new int[observations];

        for (int i = 1; i <= observations; i++) {
            O[i-1] = Integer.parseInt(oTemp[i]);
        }
        double alphaSum = alphaPass(A, B, pi, O);
        System.out.println(alphaSum);
    }

    private static double[][] matrixVectorMultiplication(double[][] vector, double[][] matrix) {

        int vectorCols = vector[0].length;
        int matrixCols = matrix[0].length;
        int matrixRows = matrix.length;

        double[][] newVector = new double[1][matrixCols];

        System.err.println("Vektorlängd: " + vectorCols + " Matriskolumner: " + matrixCols + " Rader: " + matrixRows);

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

    private static double[][] matrixInputReader(BufferedReader bReader) throws IOException {
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

    private static double[][] dotMultiply(double[][] vector1, double[][] vector2) {

        double[][] productVector = new double[1][vector1[0].length];

        if (vector1[0].length != vector2[0].length) {
            System.err.println("Cannot do dot product operation, vector lengths do not match.");
        } else {
            for (int i = 0; i < vector1[0].length; i++) {
                productVector[0][i] = vector1[0][i] * vector2[0][i];
            }
        }

        return productVector;
    }

    private static double dotProduct(double[][] vector1, double[][] vector2) {
        double[][] productVector = dotMultiply(vector1, vector2);

        double dotSum = 0;

        for (int i = 0; i < vector1[0].length; i++) {
            dotSum += productVector[0][i];
        }
        return dotSum;

    }

    private static double[][] returnColumn(double[][] matrix, int colNo) {
        double[][] vector = new double[1][matrix.length];

        for (int i = 0; i < matrix.length; i++) {
            vector[0][i] = matrix[i][colNo];
        }

        return vector;
    }


    private static double alphaPass(double[][] A, double[][] B, double[][] pi, int[] obsVector) {
        int firstObs = obsVector[0];

        double[][] bColumn = returnColumn(B, firstObs);
        // FUNKAR EJ
        System.err.println(Arrays.toString(bColumn[0]));

        double[][] alphaVector = dotMultiply(bColumn, pi);
        System.err.println(Arrays.toString(alphaVector[0]));

        for (int i = 1; i < obsVector.length; i++) {
            alphaVector = alphaPassIteration(alphaVector, A, B, obsVector[i]);
            System.err.println(Arrays.toString(alphaVector[0]));
        }

        double alphaSum = 0;

        for (int i = 0; i < pi[0].length; i++) {
            alphaSum += alphaVector[0][i];
        }
        return alphaSum;
    }

    private static double[][] alphaPassIteration(double[][] alphaPrev, double[][] A, double[][] B, int obsNo) {
        double[][] tempVector = new double[1][A[0].length];

        for (int i = 0; i < A[0].length; i++) {
            double[][] ACol = returnColumn(A, i);
            tempVector[0][i] = dotProduct(ACol, alphaPrev);
        }

        double[][] observationCol = returnColumn(B, obsNo);

        return dotMultiply(tempVector, observationCol);
    }
}
