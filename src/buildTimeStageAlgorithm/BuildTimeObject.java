package buildTimeStageAlgorithm;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import GeneticAlgorithm.Computation;
import GeneticAlgorithm.DataCenter;
import GeneticAlgorithm.Dataset;

public class BuildTimeObject {
	private int[][] dependencyMatrix = null;
	private int[][] BEAMatrix = null;
	List<Dataset> datasets = null;
	List<DataCenter> dataCenters = null;
	int[][] placementMatrix = null;
	
	public BuildTimeObject(List<Dataset> datasets, List<DataCenter> dataCenters,List<Computation> computations) {
		this.datasets = datasets;
		this.dataCenters = dataCenters;
		startAlgorithm(datasets,dataCenters,computations);
	}
	
	private void startAlgorithm(List<Dataset> datasets, List<DataCenter> dataCenters,List<Computation> computations) {
		dependencyMatrix = createDependencyMatrix(datasets,computations);
		int nrDatasets = datasets.size();
//		dependencyMatrix = new int[][] {   // diagonal 0, symetric
//			{2,1,2,0,0},
//			{1,3,1,2,1},
//			{2,1,2,0,0},
//			{0,2,0,2,1},
//			{0,1,0,1,2},
//		};
		BEAMatrix = calculateBEA(dependencyMatrix);
		displayMatrix(BEAMatrix);
		int[] datasetsIDS = new int[nrDatasets];
		for(int i=0;i<nrDatasets;i++) {
			datasetsIDS[i] = datasets.get(i).getIndex();
		}
		binaryPartition(BEAMatrix, BEAMatrix.length,datasetsIDS);
		placementMatrix = createPlacementMatrix();
		displayMatrix(placementMatrix);
		System.out.println();
	}
	
	private List<int[][]> binaryPartition(int[][] matrix,int n,int[] datasetsIds){
		List<int[][]> partitions = new ArrayList<int[][]>();
		int max = Integer.MIN_VALUE;
		int bestP = -1;
		for(int p=0;p<n;p++) {
			int pm = calculatePM(matrix, p);
			if(pm >= max) {
				max = pm;
				bestP = p;
			}
		}
		
		int[][] matrix1 = new int[bestP][n];
		int[][] matrix2 = new int[n-bestP][n];
		int[] m1Datasets = new int[bestP];
		int[] m2Datasets = new int[n-bestP];
		int count = 0;
		for(int i=0;i< bestP;i++) {
			for(int j=0;j<n;j++) {
				matrix1[i][j] = matrix[i][j];
			}
			m1Datasets[count] = datasetsIds[count];
			count++;
		}
		int k = 0;
		for(int i=bestP;i<n;i++) {			
			for(int j=0;j<n;j++) {
				matrix2[k][j] = matrix[i][j];
			}
			
			m2Datasets[k] = datasetsIds[count];
			k++;
			count++;
		}
		
		if(!canDistributeDatasets(m1Datasets)) {
			binaryPartition(matrix1, matrix1.length, m1Datasets);
		}
		
		if(!canDistributeDatasets(m2Datasets)) {
			binaryPartition(matrix2, matrix2.length, m2Datasets);
		}
		
		return partitions;
	}
	
	private boolean canDistributeDatasets(int[] datasetsIds) {
		DataCenter bestDc = null;
		long bestRemainingCapacity = Long.MIN_VALUE;
		long datasetsCapacity = 0;
		for(int i=0;i<datasetsIds.length;i++) {
			datasetsCapacity += getDatasetById(datasetsIds[i]).getFileSize();
		}
		
		for(DataCenter dc : dataCenters) {
			long availableSize = dc.getCapacity() - dc.getUsedCapacity();
			if(availableSize - datasetsCapacity > 0 && availableSize - datasetsCapacity > bestRemainingCapacity) {
				bestRemainingCapacity = availableSize - datasetsCapacity;
				bestDc = dc;
			}
		}
		
		if(bestDc != null) {
			distributeDatasetsToDatacenter(datasetsIds,bestDc);
			return true;
		}
		
		return false;
	}
	
	
	private void distributeDatasetsToDatacenter(int[] datasetsIds,DataCenter dc) {
		List<Dataset> datasetsToDistribute = new ArrayList<Dataset>();
		for(int i=0;i< datasetsIds.length;i++) {
			datasetsToDistribute.add(getDatasetById(datasetsIds[i]));
		}
		
		dc.setDatasetsToProcess(datasetsToDistribute);
		
	}
	
	private Dataset getDatasetById(int id) {
		return datasets.stream().filter(s -> s.getIndex() == id).findFirst().get();
	}
	
	private int[][] createPlacementMatrix(){
		int n = datasets.size();
		int m = dataCenters.size();
		int[][] matrix = new int[n][m];
		
		for(int i=0;i<datasets.size();i++) {
			Dataset d = datasets.get(i);
			for(int j=0;j<dataCenters.size();j++) {
				DataCenter dc = dataCenters.get(j);
				if(dc.getDatasetsToProcess().contains(d)) {
					matrix[i][j] = 1;
				}
			}
		}
		
		return matrix;
	}
	
	private int calculatePM(int[][] matrix, int p) {
		int pm = 0;
		int a1 =0;
		int a2 =0;
		int a3 = 0;
		int n = matrix.length;
		for(int i=0;i<p;i++) {
			for(int j=0; j<p;j++) {
				a1 += matrix[i][j];
			}
		}
		//poate doar de la p
		for(int i=p+1;i<n;i++) {
			for(int j=p+1;j<n;j++) {
				a2 += matrix[i][j];
			}
		}
		
		for(int i=0;i<p;i++) {
			for(int j=p+1;j<n;j++) {
				a3 += matrix[i][j];
			}
		}
		pm = a1 * a2 - (a3*a3);	
		
		return pm;
	}
	
	
	private int[][] createDependencyMatrix(List<Dataset> datasets, List<Computation> computations){
		int[][] dependencyMatrix = new int[datasets.size()][datasets.size()];
		for(int i=0;i < datasets.size();i++) {
			for(int j=i;j < datasets.size();j++) {
				int dependencyij = 0;
				for(Computation c : computations) {
					List<Integer> processedDatasets = c.getProcessedDatasets();
					if(processedDatasets.contains(datasets.get(i).getIndex()) && processedDatasets.contains(datasets.get(j).getIndex())) {
						dependencyij++;
					}
				}
				dependencyMatrix[i][j] = dependencyMatrix[j][i] = dependencyij;
				
			}
		}
		
		return dependencyMatrix;
	}
	
	private int[][] calculateBEA(int[][] dependencyMatrix){
		int n = dependencyMatrix.length;
		int[][] beaMatrix = new int[n][n];
		//add the first 2 columns from dependency matrix
		int index =0;
		for(int j=0;j< 2;j++) {
			for(int i=0;i< n;i++) {
				beaMatrix[i][j] = dependencyMatrix[i][j];
			}
		}
		List<Pair> pairs = new ArrayList<Pair>();
		index = 2;
		while(index < n) {
			int[] currentColumn  = getColumn(index, dependencyMatrix, n);
			int position = calculateAllContributions(beaMatrix, index, currentColumn);
			beaMatrix = insertColumn(beaMatrix,position,currentColumn,n);
			pairs.add(new Pair(index,position));
			index++;
		}
		//displayMatrix(beaMatrix);
		for(Pair p: pairs) {
			if(p.getIndex1()!=p.getIndex2()) {
			beaMatrix = swapRows(beaMatrix, p.getIndex1(), p.getIndex2(), n);
			}
		}
		//displayMatrix(beaMatrix);
		
 		return beaMatrix;
	}
	
	public int[] getColumn(int index,int[][] matrix,int length) {
		int[] column = new int[length];

		if(index == -1 || index >= length) {
			return column;
		}
		
		for(int i=0;i<length;i++) {
			column[i] = matrix[i][index];
		}

		return column;
	}
	
	public int[] getRow(int index,int[][] matrix,int length) {
		int[] row = new int[length];
		
		for(int i=0;i< length;i++) {
			row[i] = matrix[index][i];
		}
		
		return row;	
	}
	
	public int[][] getDependencyMatrix(){
		return dependencyMatrix;
	}
	
	public static void displayMatrix(int[][] matrix) {
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[1].length; j++) {
				System.out.print(matrix[i][j] + " ");
			}
			System.out.println();
		}
		System.out.println();
	}
	
	//calculate best contribution and return the position where the column is best placed;
	private int calculateAllContributions(int[][] matrix,int index,int[] column){
		int n = matrix.length;
		int max = Integer.MIN_VALUE;
		int position = -2;
		for(int i = -1;i< index;i++) {
			int a1[] = getColumn(i, matrix, n);
			int a2[] = column;
			int a3[] = getColumn(i+1, matrix, n);
			int result = contribution(a1,a2,a3,n);
			//DACA NU MERE DIN PRIMA CA IN EXEMPLU SCOS = 
			if(result >= max) {
				max = result;
				position = i+1;
			}
		}
		
		return position;
		
	}
	
	private int contribution(int[] a1,int a2[],int a3[],int length) {
		int contribution = 0;
		contribution = 2*bond(a1,a2,length) + 2*bond(a2,a3,length) - 2*bond(a1,a3,length);
		
		return contribution;
	}
	
	private int bond(int col1[],int col2[],int length) {
		int result = 0;
		for(int i=0;i< length; i++) {
			result+= col1[i] * col2[i];
		}
		
		return result;
	}
	
	//insert the column at the position
	private int[][] insertColumn(int[][] matrix,int position,int[] columnToInsert,int n){
		int[][] tmpMatrix = new int[n][n];
		for(int i=0;i<n;i++) {
			for(int j=0;j<n;j++) {
				tmpMatrix[i][j] = matrix[i][j];
			}
		}
		//displayMatrix(matrix);
		//until the index everything is the same
		//after index shift columns by one
		for(int j=n-1;j>position;j--) {
			for(int i=0;i<n;i++) {
				matrix[i][j] = matrix[i][j-1];
			}
		}
		
		for(int i=0;i< n;i++) {
			matrix[i][position] = columnToInsert[i];
		}

		//displayMatrix(matrix);
		return matrix;
	}
	
	private int[][] swapRows(int[][] matrix,int index1,int index2,int n){
		int[][] tmpMatrix = new int[n][n];
		for(int i=0;i<n;i++) {
			for(int j=0;j<n;j++) {
				tmpMatrix[i][j] = matrix[i][j];
			}
		}
		
		for(int i=0;i<n;i++) {
			matrix[index1][i] = tmpMatrix[index2][i];
		}
		for(int i=0;i<n;i++) {
			matrix[index2][i] = tmpMatrix[index1][i];
		}
		
		return matrix;
	}
	
	public int[][] getPlacementMatrix(){
		return this.placementMatrix;
	}
//	private int contribution(int[][] matrix,int length) {
//		int result = 0;
//		for(int i=0;i<length;i++) {
//			for(int j=0;j<length;j++) {
//				int a1 = matrix[i][j];
//				int a2 = 0;
//				int a3 = 0;
//				if(!(j-1 < 0)) {
//					a2 = matrix[i][j-1];
//				}
//				if(!(j+1 >= length)) {
//					a3 = matrix[i][j+1];
//				}
//				result += 	a1*(a2+a3);
//			}
//			
//		}
//		
//		return result;
//	}
//	
}
