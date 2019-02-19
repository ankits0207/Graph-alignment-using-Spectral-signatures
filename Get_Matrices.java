import java.util.HashMap;

public class Get_Matrices {
	
	public int[][] get_adjacency_matrix(Graph g, HashMap<String, Integer> hm){
		System.out.println("Computing adjacency matrix...");
		int[][] adjacency_matrix = new int[g.al_of_nodes.size()][g.al_of_nodes.size()];
		
		for(MyEdge me : g.al_of_edges){
			adjacency_matrix[hm.get(me.source)][hm.get(me.target)] = 1;
			adjacency_matrix[hm.get(me.target)][hm.get(me.source)] = 1;
		}		
		
		return adjacency_matrix;		
	}
	
	public int[][] get_degree_matrix(int[][] adjacency_matrix, HashMap<String, Integer> hm){
		System.out.println("Computing degree matrix...");
		HashMap<String, Integer> degree_map = new HashMap<String, Integer>();
		int[][] degree_matrix = new int[adjacency_matrix.length][adjacency_matrix.length];
		for(int i = 0; i < adjacency_matrix.length; i++)
		{
			int temp_sum = 0;
			for(int j = 0; j < adjacency_matrix.length; j++)
			{
				temp_sum += adjacency_matrix[i][j];
			}
			degree_map.put(Integer.toString(i), temp_sum);
		}
		for(String i : degree_map.keySet())
		{
			degree_matrix[Integer.parseInt(i)][Integer.parseInt(i)] = degree_map.get(i);			
		}
		return degree_matrix;
	}
	
	public int[][] get_identity_matrix(int[][] adjacency_matrix){
		System.out.println("Computing identity matrix...");
		int n = adjacency_matrix.length;
		int[][] identity_matrix = new int[n][n];
		for(int i = 0; i < n; i++)
		{
			identity_matrix[i][i] = 1;
		}
		return identity_matrix;		
	}
	
	public double[][] get_multiplied_matrices(double[][] a, double[][] b, int size){
		System.out.println("Performing matrix multiplication...");
		double[][] result_mat = new double[size][size];
		for(int i = 0; i < size; i++)
		{
			for(int j = 0; j < size; j++)
			{
				for(int k = 0; k < size; k++)
				{
					result_mat[i][j] += a[i][k]*b[k][j];
				}
			}
		}
		return result_mat;
	}
	
	/*public double[][] get_normalized_laplacian_matrix(int[][] adjacency_mat, int[][] degree_mat, int[][] identity_mat){
		System.out.println("Computing normalized Laplacian matrix...");
		int n = adjacency_mat.length;
		int m = degree_mat.length;
		double[][] i_minus_w = new double[n][n];
		double[][] modified_deg_mat = new double[m][m];

		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++)
			{
				i_minus_w[i][j] = identity_mat[i][j] - adjacency_mat[i][j];
			}
		}
		
		for(int i = 0; i < m; i++)
		{
			for(int j = 0; j < m; j++)
			{
				modified_deg_mat[i][j] = Math.sqrt(degree_mat[i][j]);
			}
		}
		
		double[][] d1 = get_multiplied_matrices(modified_deg_mat, i_minus_w, n);
		double[][] d2 = get_multiplied_matrices(d1, modified_deg_mat, n);
		return d2;		
	}*/
	
	public double[][] get_normalized_laplacian_matrix(int[][] adj_matrix, HashMap<String, Integer> deg_map){
		System.out.println("Computing normalized Laplacian matrix...");
		double[][] normalized_laplacian_matrix = new double[adj_matrix.length][adj_matrix.length];
		for(int i = 0; i < adj_matrix.length; i++)
		{
			for(int j = 0; j < adj_matrix.length; j++)
			{
				if(i == j && deg_map.get(Integer.toString(j)) != 0)
				{
					normalized_laplacian_matrix[i][j] = 1;
				}
				else if(adj_matrix[i][j] == 1)
				{
					normalized_laplacian_matrix[i][j] = -1/Math.sqrt(deg_map.get(Integer.toString(i))*deg_map.get(Integer.toString(j)));
				}
				else
				{
					normalized_laplacian_matrix[i][j] = 0;
				}				
			}
		}
		return normalized_laplacian_matrix;		
	}
}
