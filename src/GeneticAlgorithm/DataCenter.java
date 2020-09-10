package GeneticAlgorithm;
import java.util.ArrayList;
import java.util.List;

public class DataCenter {

	private int index;
	private long capacity;
	private long usedCapacity = 0;
	private double estimatedBandwidth;
	private int[] latencyVector;
	private double estimatedProcessingCapability;    // MB/sec
	private List<Dataset> datasetsToProcess;
	
	public DataCenter(int index,long capacity,double estimatedBandiwidth,double estimatedProcessingCapability,int[] latencyVector) {
		this.index = index;
		this.capacity = capacity;
		this.estimatedBandwidth = estimatedBandiwidth;
		this.latencyVector = latencyVector;
		this.estimatedProcessingCapability = estimatedProcessingCapability;
		this.datasetsToProcess = new ArrayList<Dataset>();
	}

	public int getIndex() {
		return index;
	}

	public void setIndex(int index) {
		this.index = index;
	}

	public long getCapacity() {
		return capacity/1000000;
	}

	public void setCapacity(long capacity) {
		this.capacity = capacity;
	}

	public double getEstimatedBandwidth() {
		return estimatedBandwidth;
	}

	public void setEstimatedBandwidth(double estimatedBandwidth) {
		this.estimatedBandwidth = estimatedBandwidth;
	}

	public int[] getLatencyVector() {
		return latencyVector;
	}

	public void setLatencyVector(int[] latencyVector) {
		this.latencyVector = latencyVector;
	}

	public double getEstimatedProcessingCapability() {
		return estimatedProcessingCapability * 100;    //suppose 2.2Ghz -> 100mb in 0.5 sec
	}

	public void setEstimatedProcessingCapability(double estimatedProcessingCapability) {
		this.estimatedProcessingCapability = estimatedProcessingCapability;
	}

	public List<Dataset> getDatasetsToProcess() {
		return datasetsToProcess;
	}

	public void setDatasetsToProcess(List<Dataset> datasetsToProcess) {
		this.datasetsToProcess = datasetsToProcess;
		for(Dataset d : datasetsToProcess) {
			usedCapacity += d.getFileSize();
		}
	}

	public long getUsedCapacity() {
		return usedCapacity;
	}

	public void setUsedCapacity(long usedCapacity) {
		this.usedCapacity = usedCapacity;
	}
	
	
	
}
