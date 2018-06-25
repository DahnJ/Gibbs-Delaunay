
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/point_generators_3.h>
// #include <CGAL/Kernel/global_functions.h>
#include <vector>
#include <iostream>
#include <math.h>
#include <random>
#include <time.h>
#include <tuple>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel	K;
typedef K::FT							Weight;
typedef K::Point_3						Point;
typedef K::Weighted_point_3					Weighted_point;

typedef std::tuple<Point,Point,Point> Face;
typedef std::tuple<Point,Point> 		 Edge;

typedef K::Vector_3						Vector;
typedef CGAL::Creator_uniform_3<double, Point> 			Creator;


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_cell_base_3<K>   Cbb;
typedef CGAL::Regular_triangulation_vertex_base_3<K>         Vb;
typedef CGAL::Regular_triangulation_cell_base_3<
           K, Cbb, CGAL::Discard_hidden_points>               Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>          Tds;
typedef CGAL::Regular_triangulation_3<K,Tds>                 Rt;


typedef Rt::Tetrahedron						Tetrahedron;

typedef Rt::Cell_handle						Cell_handle;


double minimumEdgeLength( const Tetrahedron& t){
	double minimum_edge_length = 1000000;
	for (int i = 0; i<3; ++i){
		for (int j = i+1; j<4; ++j){
			double edge_length = CGAL::squared_distance(t.vertex(i), t.vertex(j));
			if (edge_length < minimum_edge_length) {
				minimum_edge_length = edge_length;
			}
		}
	}
	return CGAL::sqrt(minimum_edge_length);
}

double edgeLength( const Edge & e ) {
	return(CGAL::sqrt(CGAL::squared_distance(std::get<0>(e),std::get<1>(e))));

}

double faceArea( const Face & f) {
	return(CGAL::sqrt(CGAL::squared_area(std::get<0>(f), std::get<1>(f), std::get<2>(f))));
}


double circumradius( const Tetrahedron& t) {
	return CGAL::sqrt(CGAL::squared_distance(CGAL::circumcenter(t),t.vertex(0)));
}

// K::FT volume() const { // method of Tetrahedron that returns the signed volume }

double surfaceArea( const Tetrahedron& t) {
	double surface_area = 0;
	for ( int i = 0; i<4; ++i){
		surface_area += CGAL::sqrt(CGAL::squared_area(t.vertex(i), t.vertex(i+1), t.vertex(i+2)));
	}
	return surface_area;
}



bool isWithinUnitBox(const Point& p) {
	return( (p.x() < 1) && (p.x() >= 0)
	     && (p.y() < 1) && (p.y() >= 0)
	     &&	(p.z() < 1) && (p.z() >= 0));


}


bool isWithinUnitBox(const Weighted_point& p) {
	return(isWithinUnitBox(p.point()));
}


double bounceBack(double r){
	double result;
	if (r>1) { result = 2 - r; }
	else if (r < 0) { result = -r; }
	else { result = r; }
}


Point bouncePoint(const Point& p){
	Point  mp(bounceBack(p.x()), bounceBack(p.y()), bounceBack(p.z()));
	return mp;
}



std::random_device	rand_dev;
std::mt19937		generator(rand_dev());
//std::mt19937		generator(5);
std::uniform_real_distribution<>	unif(0.0,1.0);
std::normal_distribution<> 		norm{0,0.1};


double uniformDistribution( double max = 1.0 ){
	std::uniform_real_distribution<>  uniform(0.0,max);
	return(uniform(generator));
}

Vector normalDistributionVector( double std = 0.01 ){
	return(Vector(norm(generator),norm(generator),norm(generator)));
}

Point uniformDistributionPoint(){
	return(Point(unif(generator),unif(generator),unif(generator)));
}

// TODO: Weight 
Weighted_point uniformDistributionWeightedPoint( double max = 0.01 ){
	Point p = uniformDistributionPoint();
	Weight w = uniformDistribution(max);
	return(Weighted_point(p,w));
}


// TODO: Polymorphism for Power - Delaunay?
// Class preparation
// Constructor
// 	- one more parameter
// double updateEnergy()
// double expEnergy()
// Vertex_Handle chooseRandomVertex()
// 	- different return (Dt vs Rt)
// void iNitialize()
// 	- weighted points vs points
// add, remove, move
// 	- just checking for empty handle + different handles + weights
// step
// 	- again, checking for empty hadle + different handles
// void iterate()
class Gibbs_Delaunay{
private:
	double minimum_edge_length;
	double maximum_circumradius;
	double theta;
	double intensity;
	
	double energy = 0;
	int forbidden; // number of forbidden cells/edges detected 

	int number_of_active_points = 0;

	double max_weight;

	Rt T;


public:
	Gibbs_Delaunay( double _minimum_edge_length = 0.01, double _maximum_circumradius = 0.1, double _theta = 1.0, double _intensity = 1.0, double _max_weight = 0.01):
		minimum_edge_length(_minimum_edge_length), maximum_circumradius(_maximum_circumradius), theta(_theta), intensity(_intensity), max_weight(_max_weight), forbidden(0) {  }; 

	int numberOfPoints() const {
		return T.number_of_vertices();	
	}
	

	bool isActive(const Cell_handle& c) const {
		return(isWithinUnitBox(CGAL::centroid(T.tetrahedron(c))));

	}

	// TODO: Don't have the same function twice
	double updateEnergy(const Rt::Finite_cells_iterator& begin, const Rt::Finite_cells_iterator& end, bool add = true, bool update = true) {
		double energy_update = 0;

		for (Rt::Finite_cells_iterator cell = begin; cell != end; ++cell) {
			if (T.is_infinite(cell)) { continue; }
			Tetrahedron t = T.tetrahedron(cell);
			if (!isActive(cell)) { continue; } 
			if (minimumEdgeLength(t) < minimum_edge_length){
				std::cout << "Short edge: " << minimumEdgeLength(t)  << std::endl;
				if (update) {
					if (add) {
						forbidden++;
					}
					else {
						forbidden--;
					}
				}
			}
			if (circumradius(t) > maximum_circumradius){
				std::cout << "Large circumradius: " << circumradius(t)  << std::endl;
				if (update) {
					if (add) {
						forbidden++;
					}
					else {
						forbidden--;
					}
				}
			}

			energy_update += theta*surfaceArea(t);

		}
	
		// Update energy
		if (update){
			if (add){
				energy += energy_update;
			}
			else {
				energy -= energy_update;
			}
		}

		return(energy_update);
	}
	
	double updateEnergy(const std::vector<Cell_handle>& cells, bool add = true, bool update = true) {
		double energy_update = 0;
		for (const Cell_handle& cell: cells){
			if (T.is_infinite(cell)) { continue; }
			Tetrahedron t = T.tetrahedron(cell);
			if (!isActive(cell)) { continue; } 
			if (minimumEdgeLength(t) < minimum_edge_length){
				if (update) {
					if (add) {
						forbidden++;
					}
					else {
						forbidden--;
					}
				}

			}
			if (circumradius(t) > maximum_circumradius){
				if (update) {
					if (add) {
						forbidden++;
					}
					else {
						forbidden--;
					}
				}
			}

			energy_update += theta*surfaceArea(t);

		}
		
		// Update energy
		if (update) {
			if (add){
				energy += energy_update;
			}
			else {
				energy -= energy_update;
			}
		}	

		return(energy_update);
	}


	double expEnergy(double energy){
		return(exp(-0.01*energy));
	}

	Rt::Vertex_handle chooseRandomVertexBetter(){
		// Cell_handle containing_cell = T.locate(p);
		// Tetrahedron t = T.tetrahedron(containing_cell);
		// Rt::Vertex_handle v;
		// for (int i = 0; i < 3; ++i){
		// 	if (isWithinUnitBox(t.vertex(i))) { v = t.vertex(i); break; }
		// }	
		bool done = false;
		Rt::Vertex_handle v;
		while (!done) {
			Point p = uniformDistributionPoint();
			v = T.nearest_power_vertex(p);
			if (isWithinUnitBox(T.point(v))) { done = true; }

		}
		return v;
		
			
	}
	
	// TODO: Improve this
	Rt::Vertex_handle chooseRandomVertex(){
		int n = T.number_of_vertices();
		std::uniform_int_distribution<> sample(0,n-1);
		bool sampling_done = false;
		Rt::Vertex_handle vertex;
		while (!sampling_done){
			int random_point = sample(generator);
			// std::cout << "Random number: " << random_point << std::endl;	
			// int random_point = 547;
			// std::cout << "Random vertex: " << random_point << std::endl;
			int counter = 0;
			for (Rt::Finite_vertices_iterator v = T.finite_vertices_begin(); v != T.finite_vertices_end(); ++v){
				if (counter == random_point){
					vertex = v;
					break;
				}
				counter++;
			}
			// std::cout << "Chosen point: " << T.point(vertex)  << std::endl;
			if (isWithinUnitBox(T.point(vertex))) { sampling_done = true; }
		}		
		return(vertex);
	}
		

	int numberOfActivePoints() const {
		int count = 0;	
		for (Rt::Finite_vertices_iterator v = T.finite_vertices_begin(); v != T.finite_vertices_end(); ++v){
			if (isWithinUnitBox(T.point(v))) { count++; };	
		}
		return count;
	}


	void initialize(bool from_file = true) {
		if (from_file){
			std::ifstream f("initial.txt");
			f.precision(17);
			f >> T;
			std::cout << T.is_valid() << std::endl;
		} 
		else
		{
			int number_of_points = 6000;
			std::vector<Point> points;
			std::vector<Weighted_point> weighted_points;
			points.reserve(number_of_points);

			CGAL::points_on_cube_grid_3( 1, number_of_points, std::back_inserter(points), Creator() );
			for (Point& p: points){
				p += Vector(0.5,0.5,0.5);
				Weight w = 0.005;
				weighted_points.push_back(Weighted_point(p,w));
			}

			T.insert(weighted_points.begin(), weighted_points.end());
		}

		number_of_active_points = numberOfActivePoints();

		updateEnergy(T.finite_cells_begin(), T.finite_cells_end());	
		
		std::cout << "Initialization done. Initial tessellation with " << T.number_of_vertices() << " vertices and " << T.number_of_cells() << " cells created. Initial energy: " << ( forbidden ? -1 : energy ) << std::endl;

		if (forbidden){
			std::cout << "Tesselation is not permissible! Number of forbidden elements: " << forbidden << std::endl;
		}

			
	}


	double realEnergy(){
		return(updateEnergy(T.finite_cells_begin(), T.finite_cells_end(), false, false ));
	}


 	void writeHiddenPoints() const {
		std::cout << "Hidden points: " << std::endl;
		int count = 0;
		for (Rt::Finite_cells_iterator c = T.finite_cells_begin(); c != T.finite_cells_end(); ++c){
			std::vector<Weighted_point> v;
			std::copy(c->hidden_points_begin(), c->hidden_points_end(), std::back_inserter(v));
			for (Weighted_point p: v){
				std::cout << p << std::endl;
				++count;
		  	}
		}
		std::cout << "Number of hidden points: " << count;
		

	}



	Rt::Vertex_handle add(const Weighted_point& p){
		// Find containing cell
		Cell_handle containing_cell = T.locate(p);
		
		// Find cells conflicting with p - their energy will be subtracted
		std::vector<Cell_handle> conflicting_cells;
		Rt::Facet f;
		T.find_conflicts(p , containing_cell , CGAL::Oneset_iterator<Rt::Facet>(f) , std::back_inserter(conflicting_cells));
		updateEnergy(conflicting_cells, false);
		

		// Insert new point
		Rt::Vertex_handle vertex = T.insert(p, containing_cell);

		if (vertex != Rt::Vertex_handle()) { // In case the point wasn't added, i.e. vertex is an empty handle
			// Find newly created cells = cells incident to p
			std::vector<Cell_handle> incident_cells;
			T.incident_cells(vertex, std::back_inserter(incident_cells));
			updateEnergy(incident_cells);
		}			
		return(vertex);
	}

// 	Rt::Vertex_handle addRandomPoint() {
// 		bool added = false;
// 		while (!added) {
// 			Weighted_point x = uniformDistributionWeightedPoint();
// 			// Is it identical to an existing vertex?
// 
// 			// If not, try adding
// 
// 			// Find containing cell
// 			Cell_handle containing_cell = T.locate(p);
// 			
// 			// Find cells conflicting with p - their energy will be subtracted
// 			std::vector<Cell_handle> conflicting_cells;
// 			Rt::Facet f;
// 			T.find_conflicts(p , containing_cell , CGAL::Oneset_iterator<Rt::Facet>(f) , std::back_inserter(conflicting_cells));
// 			updateEnergy(conflicting_cells, false);
// 			
// 
// 			// Insert new point
// 			Rt::Vertex_handle vertex = T.insert(p, containing_cell);
// 
// 			if (vertex != Rt::Vertex_handle()) { // In case the point wasn't added, i.e. vertex is an empty handle
// 				// Find newly created cells = cells incident to p
// 				std::vector<Cell_handle> incident_cells;
// 				T.incident_cells(vertex, std::back_inserter(incident_cells));
// 				updateEnergy(incident_cells);
// 			}			
// 
// 			std::cout << "Add    ";
// 			int n = T.number_of_vertices();
// 			Weighted_point x = uniformDistributionWeightedPoint();
// 			
// 			double energy_before = energy;
// 			Rt::Vertex_handle added_vertex = add(x);
// 			
// 			std::cout << energy_before << "->" << ( forbidden ? -1 : energy ) << " "; 
// 			
// 			if (added_vertex != Rt::Vertex_handle()) {	
// 		}
// 
// 	}

	void remove(Rt::Vertex_handle vertex){
		
		// Find cells to be deleted = cells incident to p
		std::vector<Cell_handle> incident_cells;
		T.incident_cells(vertex, std::back_inserter(incident_cells));
		updateEnergy(incident_cells, false);
		Weighted_point p = T.point(vertex);

		// Remove the point
		T.remove(vertex);
		// Find containing cell
		Cell_handle containing_cell = T.locate(p);

		// Find cells conflicting with p - their energy will be subtracted
		std::vector<Cell_handle> conflicting_cells;
		Rt::Facet f;
		T.find_conflicts(p , containing_cell , CGAL::Oneset_iterator<Rt::Facet>(f) , std::back_inserter(conflicting_cells));
		updateEnergy(conflicting_cells);
		

	}
	
	// TODO: Fix the bug where energy changes after a move
	Rt::Vertex_handle move(Rt::Vertex_handle vertex_from, Weighted_point point_to){
		remove(vertex_from);
		Rt::Vertex_handle vertex = add(point_to);
		return(vertex);
	}

	void step(){
		double a = unif(generator);
		if (a < 1.0/3.0) {
			std::cout << "B, ";
			int prev_total_n = T.number_of_vertices();
			int n = number_of_active_points;
			Weighted_point x = uniformDistributionWeightedPoint();
			
			std::cout << x << ", , ";

			double energy_before = energy;
			Rt::Vertex_handle added_vertex = add(x);
		        	
			std::cout << energy_before << ", " << ( forbidden ? -1 : energy ) << ", "; 
			
			if (added_vertex != Rt::Vertex_handle()) {	
				// See if we accept the proposal
				double b = unif(generator);
				std::cout << b << ", ";
				double ratio = ( expEnergy(energy)*intensity / ( (n+1) * expEnergy(energy_before)));
				std::cout << ratio << ", ";
				if ( forbidden || (b > ratio )){ // If not accepted, roll back
					T.remove(added_vertex);
					energy = energy_before;
					forbidden = 0;
					std::cout << "Rejected, ";
					std::cout << T.number_of_vertices() << ", ";
				}
				else{
					std::cout << "Proposal accepted" << ", ";
					std::cout << T.number_of_vertices() << ", ";
				}
			}

			number_of_active_points += (T.number_of_vertices() - prev_total_n);
		}
		else if ( a > 2.0/3.0) {
			std::cout << "D, ";	
			// Choose a random point to remove
			Rt::Vertex_handle vertex = chooseRandomVertexBetter();
			Weighted_point removed_point = T.point(vertex);
			
			std::cout << removed_point << ", , ";

			int prev_total_n = T.number_of_vertices();
			int n = number_of_active_points;
			double energy_before = energy;	
			remove(vertex);
			

			std::cout << energy_before << ", " << ( forbidden ? -1 : energy ) << ", "; 
			
			// See if we accept the pproposal
			double b = unif(generator);
			std::cout << b << ", ";
			double ratio = n*expEnergy(energy) / (expEnergy(energy_before)*intensity); 
			std::cout << ratio << ", ";
			if ( forbidden || ( b > ratio)) {
				T.insert(removed_point);
				energy = energy_before;
				forbidden = 0;
				std::cout << "Rejected, ";
				std::cout << T.number_of_vertices() << ", ";
			}	
			else{
				std::cout << "Proposal accepted"<< ", ";
				std::cout << T.number_of_vertices() << ", ";

			}

			number_of_active_points += (T.number_of_vertices() - prev_total_n);

		}
		else {
			std::cout << "M, ";	
			// Choose a random point to move
			Rt::Vertex_handle vertex = chooseRandomVertexBetter();
			Weighted_point old_point = T.point(vertex);
			// Generator a random vector by which the point will move
			Vector move_by = normalDistributionVector(); 
			Weighted_point x( bouncePoint(T.point(vertex).point() + move_by)  , T.point(vertex).weight());
			// std::cout << "Proposing to move the point " << std::endl << T.point(vertex) << " by " << std::endl <<  move_by << " to " << std::endl <<  x;
			
			std::cout << old_point << ", " << x << ", ";
			int prev_total_n = T.number_of_vertices();
			int n = number_of_active_points;

			double energy_before = energy;
			Rt::Vertex_handle new_vertex = move(vertex, x);

			std::cout << energy_before << ", " << ( forbidden ? -1 : energy ) << ", "; 

			
			// See if we accept the proposal
			if (new_vertex != Rt::Vertex_handle()){
				double b = unif(generator);
				std::cout << b << ", ";
				double ratio = expEnergy(energy) / expEnergy(energy_before);
				std::cout << ratio << ", ";
				if ( forbidden || (b > ratio )) {
					T.remove(new_vertex);
					T.insert(old_point);
					energy = energy_before;
					forbidden = 0;
					std::cout << "Rejected, ";
					std::cout << T.number_of_vertices() << ", ";
				}	
				else{
					std::cout << "Proposal accepted" << ", ";
					std::cout << T.number_of_vertices() << ", ";

				}
			}
			number_of_active_points += (T.number_of_vertices() - prev_total_n);

		}
	}


	int numberOfInactivePoints() const {
		int count = 0;	
		for (Rt::Finite_vertices_iterator v = T.finite_vertices_begin(); v != T.finite_vertices_end(); ++v){
			if (!isWithinUnitBox(T.point(v))) { count++; };	
		}
		return count;
	}

	void iterate(int number_of_iterations){
		clock_t start;
		for (int k = 1; k <= number_of_iterations; ++k){
		 	start = clock();
		 	std::cout << std::endl << k << ", " ;
		 	step();
		 	std::cout << (double)(clock()-start)/CLOCKS_PER_SEC << ";";
		 	// std::cout << " " << number_of_active_points << " " <<  numberOfInactivePoints() ;
		 	// std::cout << " Real energy: " << realEnergy(); 
		}
		
		std::cout << std::endl << "isValid test: " << T.is_valid() << std::endl;
		std::cout << "forbidden :  " << forbidden << std::endl;
	}



	void writeToFile(std::string filename) const {
		std::ofstream f(filename);
		f.precision(17);
		f << T;

	}

	
	void analyze() const {
		// Get only active cells
		std::vector<Tetrahedron> active_tetrahedra;
		std::set<Face> active_faces;
		std::set<Edge> active_edges;
		std::set<Rt::Vertex_handle> active_vertices;

		for (Rt::Finite_cells_iterator cell = T.finite_cells_begin(); cell != T.finite_cells_end(); ++cell) {
			if (isActive(cell)) { 
				Tetrahedron t = T.tetrahedron(cell);				

				active_tetrahedra.push_back(t); 
				for (int i = 0; i < 4 ; ++i) {
					active_vertices.insert(cell->vertex(i));
					active_faces.insert(Face(t.vertex(i), t.vertex(i+1), t.vertex(i+2)));
					for (int j = i + 1; j < 4; ++j){
						active_edges.insert(Edge(t.vertex(i),t.vertex(j)));	
					}
				}
				
			} 
		}
		

		std::vector<double> tetrahedra_volumes;
		std::vector<double> face_surfaces;
		std::vector<double> edge_lengths;
		std::vector<double> point_weights;
		std::vector<int> point_degrees;
		

		// Data to output
		// Distributions
		for(Tetrahedron t: active_tetrahedra) { tetrahedra_volumes.push_back(t.volume()); }
		for(Face f: active_faces) { face_surfaces.push_back(faceArea(f)); }
		for(Edge e: active_edges) { edge_lengths.push_back(edgeLength(e)); }
		for(Rt::Vertex_handle v: active_vertices) { 
			point_weights.push_back(T.point(v).weight());
			point_degrees.push_back(T.degree(v));
		}

		// Scalars
		int number_of_vertices = point_weights.size();
		int number_of_cells = active_tetrahedra.size();

		std::cout << tetrahedra_volumes.size() << " " << face_surfaces.size() << " " << edge_lengths.size() << " " << point_degrees.size() << " " << point_weights.size() << std::endl;

	}




};



int main() {
	Gibbs_Delaunay GD;
	GD.initialize(true);  // has to be false, file I/O buggy for larger triangulations?
	GD.iterate(1*pow(10,5));
	std::cout << GD.numberOfPoints() << " " << GD.numberOfActivePoints() << " " << GD.numberOfInactivePoints() <<  std::endl;
	GD.writeToFile("gibbs.txt");
	GD.analyze();
	GD.writeHiddenPoints();
}
