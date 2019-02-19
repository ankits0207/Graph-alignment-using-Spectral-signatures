import java.util.ArrayList;

public class Graph {
	ArrayList<MyNode> al_of_nodes;
	ArrayList<MyEdge> al_of_edges;
}
class MyNode{
	String Id;
	String label;
}
class MyEdge{
	String Id;
	String source;
	String target;
	double weight;	
}
