import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
//сохранение файликов с графиками
public class GraphsSaver {
    private static String Ahoii = null;
    static int countFiles1;

    public static void SaveG(float[] sigmaparab, String Ahoi) throws FileNotFoundException, IOException {
       
            countFiles1 = 0;
        File fff= new File(AbsPath.PathZero(Ahoii) + "/retail");
           File[] files1 = fff.listFiles();
           for (File file : files1) {
               if (file.isFile()) {
                   countFiles1++;  //number of files
               }
           }
           int nfile = countFiles1;
           String pathby = AbsPath.PathZero(Ahoii) + "/retail/00000.txt";
           int topolog = AbsPath.PathZero(Ahoii).length() + 8;
           String pathtodo = pathby.substring(0,topolog) + nfile + ".txt"; 
       
        //save
     
        try {File ff = new File(pathtodo);
       if (ff.createNewFile())
                   System.out.println(Ahoi + " -> " +  nfile + ".txt");
               else
                   System.out.println("File already exists" + " " + pathtodo);
           } catch (Exception e) {
               System.err.println(e);
           }   
for (int j = 0; j < sigmaparab.length; j++) {
    String strokaoun = sigmaparab[j] + "\n";  //ячейка для вывода 

    try(FileOutputStream fos = new FileOutputStream(pathtodo, true))
  {
      // перевод строки в байты
      byte[] buffer = strokaoun.getBytes();
        
      fos.write(buffer, 0, buffer.length);
  }
   }
    
}

public static int Numeric(int i){
    i ++;
    return i;
}

}