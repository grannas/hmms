import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;

public class HMM2 {
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

        HMM testHMM = new HMM(A, B, pi, O);
        testHMM.viterbiAlgorithm();
    }

    private static double[][] matrixInputReader(BufferedReader bReader) throws IOException {
        System.err.print("INPUT MATRIS: ");
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

        double[][] tempColumn = returnColumn(B, firstObs);
        System.err.println(Arrays.toString(tempColumn[0]));

        double[][] alphaVector = dotMultiply(tempColumn, pi);
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


    private static void saveColumn(double[][] matrix, double[][] vector, int position) {
        for (int i = 0; i < matrix.length; i++) {
            matrix[i][position] = vector[0][i];
        }
    }

    private static class HMM {
        private double[][] A;
        private double[][] B;
        private double[][] pi;
        private int[] O;
        private double[][] delta;
        private int[][] deltaID;

        private HMM(double[][] AIn, double[][] BIn, double[][] piIn, int[] OIn) {
            A = AIn;
            B = BIn;
            pi = piIn;
            O = OIn;

            delta = new double[A.length][O.length];
            deltaID = new int[A.length][O.length];
        }

        private void viterbiAlgorithm() {
            int firstObs = O[0];

            double[][] bColumn = returnColumn(B, firstObs);

            double[][] delta1 = dotMultiply(bColumn, pi);

            saveColumn(delta, delta1, 0);

            double[][] printD = returnColumn(delta, 0);

            for (int i = 1; i < O.length; i++) {
                viterbiIteration(i);
            }

            outputMostLikelySequence();
        }

        private void viterbiIteration(int iteration) {

            for (int i = 0; i < A.length; i++) {
                double[][] AColumn = returnColumn(A, i);
                double[][] deltaPrevious = returnColumn(delta, iteration - 1);
                double[][] candidateVector = dotMultiply(deltaPrevious, AColumn);

                for (int j = 0; j < candidateVector[0].length; j++) {
                    candidateVector[0][j] *= B[i][O[iteration]];
                }

                double maxProb = 0;
                int index = 0;
                for (int k = 0; k < candidateVector[0].length; k++) {
                    if (candidateVector[0][k] > maxProb) {
                        maxProb = candidateVector[0][k];
                        index = k;
                    }
                }
                delta[i][iteration] = maxProb;
                deltaID[i][iteration] = index;
            }

            double[][] printD = returnColumn(delta, iteration);
            System.err.println("Delta matrix, iteration " + iteration + ": " + Arrays.toString(printD[0]));
        }

        private void printDeltaID() {
            for (int row = 0; row < A.length; row++) {
                System.err.println(Arrays.toString(deltaID[row]));
            }
        }

        private void outputMostLikelySequence() {

            int lastLargestIndex = 0;
            double maxProb = 0;

            for (int i = 0; i < A.length; i++) {
                if (delta[i][O.length-1] > maxProb) {
                    lastLargestIndex = i;
                    maxProb = delta[i][O.length-1];
                }
            }

            if (O.length > 1) {
                printRecursionSequence(O.length-1, lastLargestIndex);
            }
            System.out.println(lastLargestIndex);
        }

        private void printRecursionSequence(int colIndex, int lastLargestIndex) {
            if (colIndex > 1) {
                System.err.println(lastLargestIndex);
                printRecursionSequence(colIndex-1, deltaID[lastLargestIndex][colIndex]);
            }
            System.out.print(deltaID[lastLargestIndex][colIndex]+ " ");
        }
    }
}
