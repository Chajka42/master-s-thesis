public class StErrCalc {
    public static float[] Erro(int inderx, float[][] dataSIN, float[] arraytokmishenSIN) {
        float [] absoluteErrorSIN = new float [2002];
        
//расчет сигма стандартного отклонения
        for (int i = 0; i < inderx-1; i++) {
            for (int j = 0; j < 2002; j++) {
               absoluteErrorSIN[j] += Math.pow((dataSIN[i][j] - arraytokmishenSIN[j]), 2);

          }}
          float drkSIN = (1.0f/(inderx-1.0f));
          for (int j = 0; j < 2002; j++) {
            absoluteErrorSIN[j] = (float) Math.sqrt( drkSIN * absoluteErrorSIN[j]);

        }
        return absoluteErrorSIN;

    }

    public static float[] ErroNeu(int inderx, float[][] dataSIN, float[] arraytokmishenSIN, float[] sigmaparab) {
        float [] aNErrorTokSIN = new float [2002];
        
        for (int i = 0; i < inderx-1; i++) {
            for (int j = 0; j < 2002; j++) {
               
               aNErrorTokSIN[j] += (float) Math.pow((dataSIN[i][j]*sigmaparab[j]/1.7955564f - arraytokmishenSIN[j]*sigmaparab[j]/1.7955564f), 2);
          }}
          float drkSIN = (1.0f/(inderx-1.0f));
          for (int j = 0; j < 2002; j++) {
            
            aNErrorTokSIN[j] = (float) Math.sqrt( drkSIN * aNErrorTokSIN[j]);
        }
        return aNErrorTokSIN;

    }
}