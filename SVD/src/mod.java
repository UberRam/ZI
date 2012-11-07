import Jama.Matrix;
import Jama.SingularValueDecomposition;


public class mod {

	public static void printMatrix(double[][] givenArray)
	{
		System.out.println();
		 for (int i=0;i<givenArray.length;i++)
		 {
			 for (int j=0;j<givenArray[i].length;j++)
			 {				
				 System.out.print(givenArray[i][j]+" \t");
			 }
			 System.out.println();
		 }
	}
	public static void printVector(double[][] givenArray, int column)
	{
		System.out.println();
		 for (int i=0;i<givenArray.length;i++)
		 {
			
				 System.out.print(givenArray[i][column]+"\t");
				 System.out.println();
		 }
	}
	
	//drukowanie macierzy trojkatnej - koncowej
	public static void printTriangleMatrix(double[][] mat)
	{
		System.out.println("drukowanie macierzy trojkatnej");
		for(int i=0; i < mat.length; i++)
		{
			for(int j=0; j<mat[i].length-1; j++)
			{
				if(i>j)
					System.out.print(mat[i][j] + "\t");
			}
			System.out.print("\n");
		}
	}
	
    public static double[][] Mnozenie(double[][] _M1, double[][] _M2){

    	 

        double[][] Wynik = new double[_M1.length][_M2[0].length]; 

        if(_M1[0].length != _M2.length){
System.out.println("èLE");

        }

        else{

 

            Wynik = new double[_M1.length][_M2[0].length]; 


            for(int i = 0; i < _M1.length; ++i){ 

                for(int j = 0; j <_M2[0].length; ++j){ 

                    for(int k = 0; k <_M2.length; ++k){ 

                        Wynik[i][j] += _M1[i][k] * _M2[k][j];

                    }

                }

            }

        }

        return Wynik;

    }

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		


		double[][] array = {
				 /*Human*/ 		{1.,0.,0.,1.,0.,0.,0.,0.,0.},
				 /*Interface*/ 	{1.,0.,1.,0.,0.,0.,0.,0.,0.},
				 /*Computer*/	{1.,1.,0.,0.,0.,0.,0.,0.,0.},
				 /*User*/ 		{0.,1.,1.,0.,1.,0.,0.,0.,0.},
				 /*System*/ 	{0.,1.,1.,2.,0.,0.,0.,0.,0.},
				 /*Response*/ 	{0.,1.,0.,0.,1.,0.,0.,0.,0.},
				 /*Time*/ 		{0.,1.,0.,0.,1.,0.,0.,0.,0.},
				 /*EPS*/ 		{0.,0.,1.,1.,0.,0.,0.,0.,0.},
				 /*Survey*/ 	{0.,1.,0.,0.,0.,0.,0.,0.,1.},
				 /*Trees*/ 		{0.,0.,0.,0.,0.,1.,1.,1.,0.},
				 /*Graph*/ 		{0.,0.,0.,0.,0.,0.,1.,1.,1.},
				 /*Minors*/ 	{0.,0.,0.,0.,0.,0.,0.,1.,1.}
				 //x 9
				 //12
				 }; 
		double[][] macierzRelacyjna = {
				 /*Human*/ 		{1.,0.,0.,1.,0.,0.,0.,0.,0.},
				 /*Interface*/ 	{1.,0.,1.,0.,0.,0.,0.,0.,0.},
				 /*Computer*/	{1.,1.,0.,0.,0.,0.,0.,0.,0.},
				 /*User*/ 		{0.,1.,1.,0.,1.,0.,0.,0.,0.},
				 /*System*/ 	{0.,1.,1.,2.,0.,0.,0.,0.,0.},
				 /*Response*/ 	{0.,1.,0.,0.,1.,0.,0.,0.,0.},
				 /*Time*/ 		{0.,1.,0.,0.,1.,0.,0.,0.,0.},
				 /*EPS*/ 		{0.,0.,1.,1.,0.,0.,0.,0.,0.},
				 /*Survey*/ 	{0.,1.,0.,0.,0.,0.,0.,0.,1.},
				 /*Trees*/ 		{0.,0.,0.,0.,0.,1.,1.,1.,0.},
				 /*Graph*/ 		{0.,0.,0.,0.,0.,0.,1.,1.,1.},
				 /*Minors*/ 	{0.,0.,0.,0.,0.,0.,0.,1.,1.}
				 //x 9
				 //12
				 }; 
		double[][] aPomniejszonaOSrednia = {
				 /*Human*/ 		{1.,0.,0.,1.,0.,0.,0.,0.,0.},
				 /*Interface*/ 	{1.,0.,1.,0.,0.,0.,0.,0.,0.},
				 /*Computer*/	{1.,1.,0.,0.,0.,0.,0.,0.,0.},
				 /*User*/ 		{0.,1.,1.,0.,1.,0.,0.,0.,0.},
				 /*System*/ 	{0.,1.,1.,2.,0.,0.,0.,0.,0.},
				 /*Response*/ 	{0.,1.,0.,0.,1.,0.,0.,0.,0.},
				 /*Time*/ 		{0.,1.,0.,0.,1.,0.,0.,0.,0.},
				 /*EPS*/ 		{0.,0.,1.,1.,0.,0.,0.,0.,0.},
				 /*Survey*/ 	{0.,1.,0.,0.,0.,0.,0.,0.,1.},
				 /*Trees*/ 		{0.,0.,0.,0.,0.,1.,1.,1.,0.},
				 /*Graph*/ 		{0.,0.,0.,0.,0.,0.,1.,1.,1.},
				 /*Minors*/ 	{0.,0.,0.,0.,0.,0.,0.,1.,1.}
				 //x 9
				 //12
				 }; 
				 double[] sredniaWektorowaDokumentow = new double[9];
				 double[] dlugosciWektorow = new double[9];
				 double[][] sumaIloczynow = new double[9][9];
				 double[][] dotProduct = new double[9][9];
				 System.out.println("array length[1] dokumenty="+array[1].length);
				 System.out.println("array length termy="+array.length);
				 for (int i=0;i<array[i].length;i++)
				 {
					 sredniaWektorowaDokumentow[i]=0;
				 }
				 System.out.print("Macierz A");
				 printMatrix(array);
				 
				 for (int i=0;i<array.length;i++)
				 {
					 for (int j=0;j<array[i].length;j++)
					 {
						 sredniaWektorowaDokumentow[j]+=array[i][j];
						
						// System.out.print(sredniaWektorowaDokumentow[j]+"\t");
					 }
					 //System.out.println();
				 }
				 System.out.println();
				 //System.out.println("sredniaWektorowaDokumentowSREDNIA");
				 for (int i=0;i<sredniaWektorowaDokumentow.length;i++)
				 {
					 sredniaWektorowaDokumentow[i]=sredniaWektorowaDokumentow[i]/array.length;
					 
					// System.out.print(sredniaWektorowaDokumentow[i]+"\t");
				 }
				 System.out.println();
				 //System.out.println("array length="+array.length);
				 Matrix A = new Matrix(array); 
				 SingularValueDecomposition svd = new SingularValueDecomposition(A);
				 Matrix S = svd.getS();
				 Matrix U = svd.getU();
				 //double[][] aPomniejszonaOSrednia=array;
				 //double[][] macierzRelacyjna=array;
				 System.out.println("Srednia wektorowa dokumentow");
				 for (int i=0;i<aPomniejszonaOSrednia.length;i++)
				 {
					 for (int j=0;j<array[i].length;j++)
					 {
						 aPomniejszonaOSrednia[i][j]-=sredniaWektorowaDokumentow[j];
						if (i==0){
							
							 System.out.print(sredniaWektorowaDokumentow[j]+"\t");
						}
						
					 }
					
				 }
				 System.out.print("\n\naPomniejszonaOSrednia");
				 printMatrix(aPomniejszonaOSrednia);
				 //printMatrix(array);
				 for (int i=0;i<9;i++)
				 {
					 for (int j=0;j<9;j++)
					 {
						 sumaIloczynow[i][j]=0;	 
					 }
					 dlugosciWektorow[i]=0;
				 }
				 System.out.println();
//				 for (int i=0;i<9;i++)
//				 {
//					 for (int j=0;j<9;j++)
//					 {
//						 System.out.print(sumaIloczynow[i][j]+"\t");
//					 }
//					 System.out.println();
//				 }
				 System.out.println("SUMA ILOCZYNOW");
				 for (int i=0;i<9;i++)
				 {
					 for (int j=0;j<9;j++)
					 {
						 for (int k=0;k<12;k++)
						 {
							 sumaIloczynow[i][j]+=(aPomniejszonaOSrednia[k][i]*aPomniejszonaOSrednia[k][j]);
							 
						 }
					 }
				 }
				 printMatrix(sumaIloczynow);
				 System.out.println();
				 System.out.println("DlugosciWektorow");
				 System.out.println();
				 for (int i=0;i<9;i++)
				 {
					 for (int j=0;j<12;j++)
					 {
						 dlugosciWektorow[i]+=Math.pow(aPomniejszonaOSrednia[j][i],2);
						 //System.out.print(aPomniejszonaOSrednia[j][i]+"\t");
					 }
					 //System.out.print(dlugosciWektorow[i]+"\t");
					 dlugosciWektorow[i]=Math.sqrt(dlugosciWektorow[i]);
					 System.out.print(dlugosciWektorow[i]+"\t");
				 }
				 System.out.println();
				 System.out.println();
//				 for (int i=0;i<9;i++)
//				 {
//					 for (int j=0;j<9;j++)
//					 {
//						 if (i+1<9)
//						 System.out.print(sumaIloczynow[i][j]+"\t");
//					 }
//					 System.out.println();
//				 }
				 System.out.println();
				 System.out.println("DOTPRODUCT patrzeÊ na indeksy!");
				 for (int i=0;i<9;i++)
				 {
					 for (int j=0;j<9;j++)
					 {
						 dotProduct[i][j]=sumaIloczynow[i][j]/(dlugosciWektorow[i]*dlugosciWektorow[j]);
						 System.out.print(dotProduct[i][j]+"\t");
					 }
					 System.out.println();
				 }
//				for(int i=0; i < macierzRelacyjna.length; i++)
//				{
//					for(int j=0; j<macierzRelacyjna[i].length-1; j++)
//					{
//						if(i>j)
//							System.out.print(macierzRelacyjna[i][j] + "\t");
//						for(int k=0;k<macierzRelacyjna[i].length-1;k++)
//						{
//							//System.out.println("traingular="+(macierzRelacyjna.length-1));
//							macierzRelacyjna[i][j]=0;
//						}
//						for(int k=0;k<macierzRelacyjna[i].length-1;k++)
//						{
//							macierzRelacyjna[i][j]=macierzRelacyjna[i][j]+aPomniejszonaOSrednia[i][k]*aPomniejszonaOSrednia[i][k+1];
//						}
//					}
//					System.out.print("\n");
//				}
//				 printTriangleMatrix(macierzRelacyjna);
				 //printVector(array,0);
				 printTriangleMatrix(dotProduct);
				 System.out.println("PO SVD");
				 Matrix V= svd.getV().transpose();
				 double [][] uKolumn = new double[12][2];
				 for (int i=0;i<12;i++)
				 {
					 uKolumn[i][0]=U.getArrayCopy()[i][0];
					 uKolumn[i][1]=U.getArrayCopy()[i][1];
				 }
				 double [][] vRow = new double[2][9];
				 for (int i=0;i<9;i++)
				 {
					 vRow[0][i]=V.getArrayCopy()[0][i];
					 vRow[1][i]=V.getArrayCopy()[1][i];
				 }
				 double [][] sMatrix = new double[2][2];
				 sMatrix[0][1]=S.getArrayCopy()[0][1];
				 sMatrix[1][1]=S.getArrayCopy()[1][1];
				 sMatrix[1][0]=S.getArrayCopy()[1][0];
				 sMatrix[0][0]=S.getArrayCopy()[0][0];
				 //System.out.println("MATRIX A");
				 //A.print(0, 2);
				 System.out.println("MATRIX V");
				 V.print(0, 2);
				 printMatrix(vRow);
				 System.out.println("MATRIX U");
				 U.print(0, 2);
				 printMatrix(uKolumn);
				 System.out.println("MATRIX S - diagonalna");
				 S.print(0, 2);
				 printMatrix(sMatrix);
				 System.out.println("____");
				 System.out.println(uKolumn.length+"sMatrix"+sMatrix[0].length);
				 double[][] macierzUS = new double[uKolumn.length][sMatrix[0].length];
				 macierzUS=Mnozenie(uKolumn, sMatrix);
				 System.out.println(macierzUS.length+"macierzUS"+macierzUS[0].length);
				 printMatrix(macierzUS);
				 double[][] macierzUSV = new double[macierzUS.length][vRow[0].length];
				 macierzUSV=Mnozenie(macierzUS, vRow);
				 printMatrix(macierzUSV);
				 System.out.println("");
				 
				 for (int i=0;i<array[i].length;i++)
				 {
					 sredniaWektorowaDokumentow[i]=0;
				 }				 
				 for (int i=0;i<array.length;i++)
				 {
					 for (int j=0;j<array[i].length;j++)
					 {
						 sredniaWektorowaDokumentow[j]+=array[i][j];
						
						// System.out.print(sredniaWektorowaDokumentow[j]+"\t");
					 }
					 //System.out.println();
				 }
				 System.out.println();
				 //System.out.println("sredniaWektorowaDokumentowSREDNIA");
				 for (int i=0;i<sredniaWektorowaDokumentow.length;i++)
				 {
					 sredniaWektorowaDokumentow[i]=sredniaWektorowaDokumentow[i]/macierzUSV.length;
					 
					// System.out.print(sredniaWektorowaDokumentow[i]+"\t");
				 }
				 System.out.println();
				 aPomniejszonaOSrednia=macierzUSV; //PRZYPISANIE!!!
//				 System.out.println("Srednia wektorowa dokumentow");
//				 for (int i=0;i<aPomniejszonaOSrednia.length;i++)
//				 {
//					 for (int j=0;j<array[i].length;j++)
//					 {
//						 aPomniejszonaOSrednia[i][j]-=sredniaWektorowaDokumentow[j];
//						if (i==0){
//							
//							 System.out.print(sredniaWektorowaDokumentow[j]+"\t");
//						}
//						
//					 }
//					
//				 }
//				 System.out.print("\n\naPomniejszonaOSrednia");
//				 printMatrix(aPomniejszonaOSrednia);
				 //printMatrix(array);
				 for (int i=0;i<9;i++)
				 {
					 for (int j=0;j<9;j++)
					 {
						 sumaIloczynow[i][j]=0;	 
					 }
					 dlugosciWektorow[i]=0;
				 }
				 System.out.println();
				 System.out.println("SUMA ILOCZYNOW");
				 for (int i=0;i<9;i++)
				 {
					 for (int j=0;j<9;j++)
					 {
						 for (int k=0;k<12;k++)
						 {
							 sumaIloczynow[i][j]+=(aPomniejszonaOSrednia[k][i]*aPomniejszonaOSrednia[k][j]);
							 
						 }
					 }
				 }
				 printMatrix(sumaIloczynow);
				 System.out.println();
				 System.out.println("DlugosciWektorow");
				 System.out.println();
				 for (int i=0;i<9;i++)
				 {
					 for (int j=0;j<12;j++)
					 {
						 dlugosciWektorow[i]+=Math.pow(aPomniejszonaOSrednia[j][i],2);
						 //System.out.print(aPomniejszonaOSrednia[j][i]+"\t");
					 }
					 //System.out.print(dlugosciWektorow[i]+"\t");
					 dlugosciWektorow[i]=Math.sqrt(dlugosciWektorow[i]);
					 System.out.print(dlugosciWektorow[i]+"\t");
				 }
				 System.out.println();
				 System.out.println();

				 System.out.println();
				 System.out.println("DOTPRODUCT patrzeÊ na indeksy!");
				 for (int i=0;i<9;i++)
				 {
					 for (int j=0;j<9;j++)
					 {
						 dotProduct[i][j]=sumaIloczynow[i][j]/(dlugosciWektorow[i]*dlugosciWektorow[j]);
						 System.out.print(dotProduct[i][j]+"\t");
					 }
					 System.out.println();
				 }
				 printTriangleMatrix(dotProduct);
	}

}
