package GeneticAlgorithm;
import java.util.ArrayList;
import java.util.List;

public class Request {

	private List<Computation> computations;
	private int[][] associationMatrix;
	private int numberOfComputations;
	
	public Request(int numberOfComputations) {
		computations = new ArrayList<Computation>();
		this.numberOfComputations = numberOfComputations;
		this.associationMatrix = new int[numberOfComputations][Utils.datasetsCount];
	}
	
	public void setAssociationMatrix(int[][] associationMatrix) {
		this.associationMatrix = associationMatrix;	
	}
	
	public  int[][] getAssociationMatrix(){
		return associationMatrix;
	}
	
	public List<Computation> getComputations(){
		return computations;
	}
	
	public void generateComputations() {
		for(int i=0;i< numberOfComputations;i++) {	
			computations.add(new Computation(i));
		}
	}
	
	
	public void createAssociationMatrix() {
		int j=0;
		for(Computation c : computations) {
			for(Integer i: c.getProcessedDatasets()) {
				associationMatrix[j][i.intValue()]= 1;
			}
			j++;
		}
	}
	
}
