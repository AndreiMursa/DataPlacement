package GeneticAlgorithm;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Computation {

	private List<Integer> processedDatasets;
	private int index;
	
	public Computation(int index) {
		this.index = index;
		processedDatasets = new ArrayList<Integer>();
		generateNeededDatasets();
	}
	
	
	private void generateNeededDatasets() {
		Random randNum = new Random();
		int nrOfNeededDatasets = 1 + randNum.nextInt(Utils.datasetsCount*Utils.datasetsPercentage/100);
		processedDatasets = new ArrayList<Integer>();
		    while (processedDatasets.size() < nrOfNeededDatasets) {
			    int nr = randNum.nextInt(Utils.datasetsCount);
		    	if(!processedDatasets.contains(nr)) {
		    		processedDatasets.add(nr);
		    	}      
		    }
	}
	
	public void printDatasetsToProcess() {
		for(Integer i:processedDatasets) {
			System.out.println(i + " ");
		}
	}
	
	
	public List<Integer> getProcessedDatasets(){
		return processedDatasets;
	}
	
	public int getIndex() {
		return index;
	}
}
