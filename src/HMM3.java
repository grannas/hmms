import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;

public class HMM3 {
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

        testHMM.BaumWelchAlgorithm();

    }

    private static void printMatrix(double[][] matrix) {
        System.err.println(matrix.length + " x " + matrix[0].length + " matrix:");
        for (int i = 0; i < matrix.length; i++) {
            System.err.println(Arrays.toString(matrix[i]));
        }
    }

    private static double[][] matrixInputReader(BufferedReader bReader) throws IOException {
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

    private static void saveColumn(double[][] matrix, double[][] vector, int position) {
        for (int i = 0; i < matrix.length; i++) {
            matrix[i][position] = vector[0][i];
        }
    }

    private static double columnSum(double[] column) {
        double sum = 0;
        for (int i = 0; i < column.length; i++) {
            sum += column[i];
        }
        return sum;
    }

    private static class HMM {
        private double[][] A;
        private double[][] B;
        private double[][] pi;
        private int[] O;

        private int T;
        private int N;
        private int M;

        private double[] alphaScaling;
        private double[] betaScaling;

        private double[][] delta;
        private int[][] deltaID;

        private double[][] alphaHistory;
        private double[][] betaHistory;

        private double[][][] diGamma;
        private double[][] altGamma;

        double logProb;

        private HMM(double[][] AIn, double[][] BIn, double[][] piIn, int[] OIn) {
            A = AIn;
            B = BIn;
            pi = piIn;
            O = OIn;

            T = O.length;
            N = B.length;
            M = B[0].length;

            alphaScaling = new double[T];
            betaScaling = new double[T];

            delta = new double[N][T];
            deltaID = new int[N][T];

            alphaHistory = new double[N][T];
            betaHistory = new double[N][T];

            diGamma = new double[N][N][T];
            altGamma = new double[N][T];

            logProb = Double.NEGATIVE_INFINITY;
        }

        private void viterbiAlgorithm() {
            int firstObs = O[0];

            double[][] bColumn = returnColumn(B, firstObs);

            double[][] delta1 = dotMultiply(bColumn, pi);

            saveColumn(delta, delta1, 0);

            double[][] printD = returnColumn(delta, 0);
            System.err.println("Delta matrisen: " + Arrays.toString(printD[0]));

            for (int i = 1; i < T; i++) {
                viterbiIteration(i);
            }

            outputMostLikelySequence();
        }

        private void viterbiIteration(int iteration) {

            for (int i = 0; i < N; i++) {
                double[][] AColumn = returnColumn(A, i);
                double[][] deltaPrevious = returnColumn(delta, iteration - 1);
                double[][] candidateVector = dotMultiply(deltaPrevious, AColumn);

                for (int j = 0; j < N; j++) {
                    candidateVector[0][j] *= B[i][O[iteration]];
                }

                double maxProb = 0;
                int index = 0;
                for (int k = 0; k < N; k++) {
                    if (candidateVector[0][k] > maxProb) {
                        maxProb = candidateVector[0][k];
                        index = k;
                    }
                }
                delta[i][iteration] = maxProb;
                deltaID[i][iteration] = index;
            }

            double[][] printD = returnColumn(delta, iteration);
            System.err.println("Delta matrisen, iteration " + iteration + ": " + Arrays.toString(printD[0]));
        }

        private void outputMostLikelySequence() {

            int lastLargestIndex = 0;
            double maxProb = 0;

            for (int i = 0; i < N; i++) {
                if (delta[i][T-1] > maxProb) {
                    lastLargestIndex = i;
                    maxProb = delta[i][T-1];
                }
            }
            if (T > 1) {
                printRecursionSequence(T-1, lastLargestIndex);
            }
            System.out.println(lastLargestIndex);
        }

        private void printRecursionSequence(int colIndex, int lastLargestIndex) {
            if (colIndex > 1) {
                printRecursionSequence(colIndex-1, deltaID[lastLargestIndex][colIndex]);
            }
            System.out.print(deltaID[lastLargestIndex][colIndex]+ " ");
        }

        private void betaPass() {
            double[][] beta = new double[1][N];
            Arrays.fill(beta[0], alphaScaling[T-1]);
            saveColumn(betaHistory, beta, T-1);

            double[][] tempBColumn;
            double[][] tempVector;
            double[][] tempA = new double[1][N];


            for (int t = T - 2; t >= 0; t--) {
                tempBColumn = returnColumn(B, O[t + 1]);

                for (int i = 0; i < M; i++) {

                    tempA[0] = A[i];
                    tempVector = dotMultiply(tempA, tempBColumn);

                    beta[0][i] = dotProduct(tempVector, returnColumn(betaHistory, t+1));
                }
                beta =  standardizeArray(beta, t);
                saveColumn(betaHistory, beta, t);
            }
        }

        private void alphaPass() {
            int firstObs = O[0];

            double[][] tempColumn = returnColumn(B, firstObs);
            double[][] alphaVector = dotMultiply(tempColumn, pi);

            alphaVector = standardizeAlpha(alphaVector, 0);

            saveColumn(alphaHistory, alphaVector, 0);

            for (int t = 1; t < T; t++) {
                alphaVector = alphaPassIteration(t);
                alphaVector = standardizeAlpha(alphaVector, t);
                saveColumn(alphaHistory, alphaVector, t);
            }
        }

        private double[][] alphaPassIteration(int time) {
            double[][] tempVector = new double[1][N];

            for (int i = 0; i < N; i++) {
                double[][] ACol = returnColumn(A, i);
                // System.err.println("ACol, iteration " + i + " : " + Arrays.toString(ACol[0]));
                // System.err.println("Alpha previous, iteration " + i + " : " + Arrays.toString(returnColumn(alphaHistory, time-1)[0]));
                tempVector[0][i] = dotProduct(ACol, returnColumn(alphaHistory, time - 1));
                //System.err.println("ACol: " + Arrays.toString(ACol[0]));
                //System.err.println("Alpha Previous: " + Arrays.toString(alphaPrev[0]));
            }

            double[][] observationCol = returnColumn(B, O[time]);

            //System.err.println("TempVect: " + Arrays.toString(tempVector[0]));
            //System.err.println("ObsCol: " + Arrays.toString(observationCol[0]));

            return dotMultiply(tempVector, observationCol);
        }

        private void updateGamma() {
            for (int t = 0; t < T; t++) {
                for (int i = 0; i < N; i++) {
                    // gamma[i][t] = alphaHistory[i][t] * betaHistory[i][t];
                }
            }
        }

        private void updateDiGamma() {
            double diGammaValue;
            double gamma;
            for (int t = 0; t < T - 1; t++) {
                for (int i = 0; i < N; i++) {
                    gamma = 0;
                    for (int j = 0; j < N; j++) {
                        diGammaValue = alphaHistory[i][t] * A[i][j] * B[j][O[t + 1]] * betaHistory[j][t + 1] / columnSum(returnColumn(alphaHistory, T -1)[0]);
                        diGamma[i][j][t] = diGammaValue;
                        gamma += diGammaValue;
                    }
                    altGamma[i][t] = gamma;
                }
            }

            double[][] lastAltGamma = returnColumn(alphaHistory, T - 1);
            double tempSum = columnSum(lastAltGamma[0]);
            for (int i = 0; i < N; i++) {
                lastAltGamma[0][i] /= tempSum;
            }
            saveColumn(altGamma, lastAltGamma, T-1);
        }

        private void reestimateParameters() {
                reestimateA();
                reestimateB();
                reestimatePi();
        }

        private void reestimateA() {
            double diGammaSum;
            double gammaSum;
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    diGammaSum = 0;
                    gammaSum = 0;

                    for (int t = 0; t < T - 1; t++) {
                        diGammaSum += diGamma[i][j][t];
                        gammaSum += altGamma[i][t];
                    }
                    A[i][j] = diGammaSum / gammaSum;
                }
            }
        }

        private void reestimateB() {
            double gammaSum;
            double gammaSumEmission;
            for (int j = 0; j < N; j++) {
                for (int k = 0; k < M; k++) {
                    gammaSum = 0;
                    gammaSumEmission = 0;
                    for (int t = 0; t < T; t++) {
                        gammaSum += altGamma[j][t];
                        if (O[t] == k) {
                            gammaSumEmission += altGamma[j][t];
                        }
                    }
                    B[j][k] = gammaSumEmission / gammaSum;
                }
            }
        }

        private void reestimatePi() {
            for (int i = 0; i < N; i++) {
                pi[0][i] = altGamma[i][0];
            }
        }

        private double logProb() {
            double logProb = 0;
            for (int t = 0; t < T; t++) {
                logProb += Math.log(1 / alphaScaling[t]);
            }
            logProb = -logProb;
            return logProb;
        }

        private void BaumWelchAlgorithm() {
            int iteration = 0;
            int maxIterations = 25;
            double newLogProb = Double.NEGATIVE_INFINITY;

            while (iteration < maxIterations && newLogProb >= logProb) {
                logProb = newLogProb;
                alphaPass();
                betaPass();
                updateGamma();
                updateDiGamma();

                reestimateParameters();

                newLogProb = logProb();

                iteration++;
            }
            outputBaumWelchAnswer();
        }

        private void outputBaumWelchAnswer() {
            printMatrixLong(A);
            printMatrixLong(B);
        }

        private void printMatrixLong(double[][] matrix) {
            int iLength = matrix.length;
            int jLength = matrix[0].length;
            System.out.print(iLength + " " + jLength);

            for (int i = 0; i < iLength; i++) {
                for (int j = 0; j < jLength; j++) {
                    System.out.print(" " + matrix[i][j]);
                }
            }
            System.out.println();
        }

        private double[][] standardizeAlpha(double[][] array, int time) {
            double scale = columnSum(array[0]);
            alphaScaling[time] = scale;
            for (int i = 0; i < array[0].length; i++) {
                array[0][i] /= scale;
            }
            return array;
        }

        private double[][] standardizeBeta(double[][] array, int time) {
            double scale = columnSum(array[0]);
            betaScaling[time] = scale;
            for (int i = 0; i < array[0].length; i++) {
                array[0][i] /= scale;
            }
            return array;
        }

        private double[][] standardizeArray(double[][] array, int time) {
            double scale = alphaScaling[time]; // ;
            for (int i = 0; i < array[0].length; i++) {
                array[0][i] /= scale;
            }
            return array;
        }
    }
}
