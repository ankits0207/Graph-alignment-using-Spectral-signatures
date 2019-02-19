import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class Utility {
	public HashMap<String, Integer> get_HashMap(Graph g){
		int i = 0;
		HashMap<String, Integer> hm = new HashMap<String, Integer>();
		for(MyNode mn : g.al_of_nodes){
			hm.put(mn.Id, i);
			i++;			
		}
		return hm;		
	}
	
	public HashMap<String, Integer> get_degree_map(int[][] adjacency_matrix, HashMap<String, Integer> hm){
		HashMap<String, Integer> degree_map = new HashMap<String, Integer>();
		for(int i = 0; i < adjacency_matrix.length; i++)
		{
			int temp_sum = 0;
			for(int j = 0; j < adjacency_matrix.length; j++)
			{
				temp_sum += adjacency_matrix[i][j];
			}
			degree_map.put(Integer.toString(i), temp_sum);
		}
		return degree_map;
	}
	
	public Graph get_1_hop_neighbourhood(Graph g_in, MyNode seed){
		Graph g_out = new Graph();
		g_out.al_of_nodes = new ArrayList<MyNode>();
		g_out.al_of_edges = new ArrayList<MyEdge>();
		g_out.al_of_nodes.add(seed);
		for(MyEdge me : g_in.al_of_edges)
		{
			if(me.source.equals(seed.Id))
			{
				g_out.al_of_nodes.add(get_node_obj(g_in, me.target));
				g_out.al_of_edges.add(me);
			}
			else if(me.target.equals(seed.Id))
			{
				g_out.al_of_nodes.add(get_node_obj(g_in, me.source));
				g_out.al_of_edges.add(me);
			}
		}
		return g_out;
	}
	
	public Graph get_k_hop_neighbourhood(Graph g, MyNode seed, int k){
		//System.out.println("Generating " + k + " hop neighbourhood...");
		Graph g_out_k_hop = new Graph();
		g_out_k_hop.al_of_nodes = new ArrayList<MyNode>();
		g_out_k_hop.al_of_edges = new ArrayList<MyEdge>();
		
		ArrayList<Graph> al_of_graph = new ArrayList<Graph>();
		ArrayList<MyNode> al_of_nodes = new ArrayList<MyNode>();
		Graph returned_graph = null;
		al_of_nodes.add(seed);
		
		HashSet<MyNode> hs_node_1 = new HashSet<MyNode>();		
		
		for(int hop_counter = 1; hop_counter <= k; hop_counter++){
			for(MyNode mn1 : al_of_nodes)
			{
				returned_graph = get_1_hop_neighbourhood(g, mn1);
				al_of_graph.add(returned_graph);
			}
			for(Graph graph : al_of_graph)
			{
				al_of_nodes.addAll(graph.al_of_nodes);
			}
			hs_node_1.addAll(al_of_nodes);
			al_of_nodes.clear();
			al_of_nodes.addAll(hs_node_1);			
		}
		for(Graph graph : al_of_graph)
		{
			for(MyNode mn : graph.al_of_nodes)
			{
				g_out_k_hop.al_of_nodes.add(mn);
			}
			for(MyEdge me : graph.al_of_edges)
			{
				g_out_k_hop.al_of_edges.add(me);
			}
		}
		HashSet<MyNode> hs_node2 = new HashSet<MyNode>();
		HashSet<MyEdge> hs_edge = new HashSet<MyEdge>();
		
		hs_node2.addAll(g_out_k_hop.al_of_nodes);
		hs_edge.addAll(g_out_k_hop.al_of_edges);
		
		g_out_k_hop.al_of_nodes.clear();
		g_out_k_hop.al_of_edges.clear();
		
		g_out_k_hop.al_of_nodes.addAll(hs_node2);
		g_out_k_hop.al_of_edges.addAll(hs_edge);
		
		return g_out_k_hop;
	}
	
	/*public Graph get_k_hop_neighbourhood(Graph g, MyNode start_vertex, int k){
		Graph gout = new Graph();
		
		//System.out.println(k);
		
		UndirectedGraph<MyNode, MyEdge> ug = new UndirectedSparseGraph<MyNode, MyEdge>();
		for(MyEdge e : g.al_of_edges)
		{
			ug.addEdge(e, get_node_obj(g, e.source), get_node_obj(g, e.target), EdgeType.UNDIRECTED);
		}
		
		for(MyNode n : g.al_of_nodes)
		{
			ug.addVertex(n);
		}
		
		//System.out.println("Test");
		
		Filter<MyNode, MyEdge> filter = new KNeighborhoodFilter<MyNode, MyEdge>(start_vertex, k, edu.uci.ics.jung.algorithms.filters.KNeighborhoodFilter.EdgeType.IN_OUT);
		edu.uci.ics.jung.graph.Graph<MyNode, MyEdge> sub_graph = filter.transform(ug);
		
		gout.al_of_nodes = new ArrayList<MyNode>();
		gout.al_of_edges = new ArrayList<MyEdge>();
		
		for(MyNode n : sub_graph.getVertices())
		{
			gout.al_of_nodes.add(n);
		}
		
		for(MyEdge e : sub_graph.getEdges())
		{
			gout.al_of_edges.add(e);
		}
		
		return gout;
	}*/
	
	public MyNode get_node_obj(Graph g, String s){
		for(MyNode mn : g.al_of_nodes)
		{
			if(mn.Id.equals(s))
				return mn;
		}
		return null;
	}
	
	public void write_score_file(ArrayList<Score> al_of_s){
		BufferedWriter bw = null;
		FileWriter fw = null;
		try {
			System.out.println("Writing to file...");
			fw = new FileWriter("Final_Scores.txt");
			bw = new BufferedWriter(fw);
			for(Score s : al_of_s)
			{
				bw.write(s.n1.Id + " " + s.n2.Id + " " + s.score);
				bw.newLine();
			}
			System.out.println("File written successfully");
		}
		catch (IOException e) {
			e.printStackTrace();
		} 
		finally {
			try {
				if (bw != null)
				{
					bw.close();					
				}
				if (fw != null)
				{
					fw.close();
				}
			}
			catch(IOException ex)
			{
				ex.printStackTrace();
			}
		}		
	}
	
	public void write_alignment_file(ArrayList<Alignment> al_of_a){
		BufferedWriter bw = null;
		FileWriter fw = null;
		try {
			System.out.println("Writing to file...");
			fw = new FileWriter("Final_Alignment.txt");
			bw = new BufferedWriter(fw);
			for(Alignment a : al_of_a)
			{
				bw.write(a.node1 + " " + a.node2);
				bw.newLine();
			}
			System.out.println("File written successfully");
		}
		catch (IOException e) {
			e.printStackTrace();
		} 
		finally {
			try {
				if (bw != null)
				{
					bw.close();					
				}
				if (fw != null)
				{
					fw.close();
				}
			}
			catch(IOException ex)
			{
				ex.printStackTrace();
			}
		}		
	}
	
	public Seed get_first_seed(ArrayList<Score> al_of_s)
	{
		Seed seed = new Seed();
		seed.score = al_of_s.get(0).score;
		for(Score s : al_of_s)
		{
			if(s.score < seed.score)
			{
				seed.score = s.score;
				seed.n1 = s.n1;
				seed.n2 = s.n2;
			}
		}
		return seed;
	}
	
}
