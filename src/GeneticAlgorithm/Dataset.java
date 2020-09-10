package GeneticAlgorithm;
import java.io.File;

public class Dataset {

	private File dataset;
	private int index;
	//the datacenter where it is currently placed
	private DataCenter datacenter;
	
	public Dataset(File dataset,int index) {
		this.dataset = dataset;
		this.index = index;
	}
	
	
	public long getFileSize() {
		return dataset.length()/1000000;
	}
	
	public int getIndex() {
		return this.index;
	}
	
	public void setIndex(int index) {
		this.index = index; 
	}
	
	public File getFile() {
		return this.dataset;
	}
	
	
	public void setDataCenter(DataCenter datacenter) {
		this.datacenter = datacenter;
	}
	
	public DataCenter getDataCenter() {
		return this.datacenter;
	}
	
}
