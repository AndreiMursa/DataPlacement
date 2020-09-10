package GeneticAlgorithm;

public class Utils {

	public static final int population =100;
	public static final int maxGeneration = 200;
	public static final double crossoverRate = 0.5;
	public static final double mutationRate = 0.2;
	public static final int replicationFactor = 1;
	public static final int dataCentersCount = 10;   //4  //15
	public static final int datasetsCount = 30;    //10  //23
	public static final int linesPerFile = 1000000; //how many lines each file should have (maximum)
	public static final String initialFilePath = "";
	
	
	//Computation
	public static final int numberOfComputations = 500;//700  //5
	public static final int datasetsPercentage = 70;   //70
	
	//Data Center
	public static final long dataCenterCapacity = 5000*1000000L;//10737418240L;//  //10GB
//	public static final int[][] latencyMatrix = new int[][] {   // diagonal 0, symetric
//		{0,40,60,80},
//		{40,0,50,20},
//		{60,50,0,55},
//		{80,20,55,0},
//	};
	public static int[][] latencyMatrix = new int[][] {   // diagonal 0, symetric
		{0,40,60,80,10,90,20,50,60,10,10,30,90,120,40},
		{40,0,50,20,10,30,90,70,30,80,80,20,10,60,40},
		{60,50,0,55,60,20,10,40,35,90,10,10,65,75,30},
		{80,20,55,0,20,10,10,80,85,20,30,60,55,30,45}, 
		{10,10,60,20,0,30,20,65,15,75,70,10,30,25,40},
		{90,30,20,10,30,0,40,45,25,30,80,95,10,20,45},
		{20,90,10,10,20,40,0,10,15,40,55,10,20,90,75},
		{50,70,40,80,65,45,10,0,15,20,75,45,30,90,65},
		{60,30,35,85,15,25,15,15,0,30,65,55,70,70,90},
		{10,80,90,20,75,30,40,20,30,0,20,10,35,95,65},
		{10,80,10,30,70,80,55,75,65,20,0,40,55,40,10},
		{30,20,10,60,10,95,10,45,55,10,40,0,65,10,25},
		{90,10,65,55,30,10,20,30,70,35,55,65,0,20,45},//
		{120,60,75,30,25,20,90,90,70,95,40,10,20,0,60},
		{40,40,30,0,45,40,45,75,65,90,65,10,25,45,0},
	};

}
