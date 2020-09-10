package GeneticAlgorithm;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.Collectors;

public class Chromosome {
	
	private int[][] placementMatrix;
	private int datasetsCount;
	private int datacentersCount;
	private double fitnessValue;
	private int totalNrDataScheduling;
	private double totalEPT;
	
	
	public Chromosome(int datasetsCount,int datacentersCount) {
		placementMatrix = new int[datasetsCount][datacentersCount];
		this.datacentersCount = datacentersCount;
		this.datasetsCount = datasetsCount;
	}
	
	public void generateInitialPlacement(int nrDatasets,int nrDataCenters,int replicaFactor) {
		for(int i = 0;i<nrDatasets;i++) {
			Random randNum = new Random();
		    List<Integer> genNrs = new ArrayList<Integer>();
		    while (genNrs.size() < replicaFactor) {
			    int nr = randNum.nextInt(nrDataCenters);
		    	if(!genNrs.contains(nr)) {
		    		genNrs.add(nr);
		    		placementMatrix[i][nr] = 1;
		    	}      
		    }
		    
		}
	}
	
	public void displayIndividual() {
		for(int i = 0;i<datasetsCount;i++) {
			for(int j=0;j<datacentersCount;j++) {
				System.out.print(placementMatrix[i][j] + " ");
			}
			System.out.println();
		}
	}
	
	
	
	public int getTotalNrDataScheduling() {
		return totalNrDataScheduling;
	}

	public void setTotalNrDataScheduling(int totalNrDataScheduling) {
		this.totalNrDataScheduling = totalNrDataScheduling;
	}

	public double getTotalEPT() {
		return totalEPT;
	}

	public void setTotalEPT(double totalEPT) {
		this.totalEPT = totalEPT;
	}

	public void setPlacementMatrix(int[][] placementMatrix) {
		this.placementMatrix = placementMatrix;
	}
	
	public int[][] getPlacementMatrix(){
		return this.placementMatrix;
	}
	
	public double getFitnessValue() {
		return fitnessValue;
	}
	
	public void setFitnessValue(double fitnessValue) {
		this.fitnessValue = fitnessValue;
	}
}
