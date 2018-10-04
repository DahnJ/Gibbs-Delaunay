
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/point_generators_3.h>
// #include <CGAL/Kernel/global_functions.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <ctime>
#include <tuple>
#include <fstream>
#include <chrono>
#include <limits>
#include <numeric> // for vector sum
#include <string>
#include <algorithm> // for min


bool COUT = false; 
bool FOUT = false;
bool DELAUNAY = false; // Forces all proposed points to have equal weight, causing the tessellation to be Delaunay
bool ANALYZE = false;


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


// TODO: Rewrite forbidden to infinity

// Abstractions
// - easily specifiable energy function and calculation
// - delaunay x regular


template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T> v){
    if ( !v.empty() ){
        out << '[';
        std::copy( v.begin(), v.end(), std::ostream_iterator<T>(out, ", ") );
        out << "]";
    }
    return out;
}



double sign(double x){
    return ( (x>0) - (x<0) );
}



///// Functions for calculating basic stats

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

double  minimumFaceArea( const Tetrahedron&  t ){
    double minimum_face_area = 1000000;
    for ( int i = 0; i < 4; ++i ) {
        double face_area_squared = CGAL::squared_area(t.vertex(i), t.vertex(i+1), t.vertex(i+2));
        if (face_area_squared < minimum_face_area) {
            minimum_face_area = face_area_squared;
        }
    } 
    return CGAL::sqrt(minimum_face_area);
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



/////////// Functions to check and remain within the unit box

bool isWithinUnitBox(const Point& p, float unit = 1) {
    float pad = (1-unit)/2;
	return( (p.x() < unit + pad) && (p.x() >= pad)
	     && (p.y() < unit + pad) && (p.y() >= pad)
	     &&	(p.z() < unit + pad) && (p.z() >= pad));


}


bool isWithinUnitBox(const Weighted_point& p, float unit = 1) {
	return(isWithinUnitBox(p.point(), unit));
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




//////////// Sampling

std::random_device	rand_dev;
std::mt19937		generator(rand_dev());
//std::mt19937		generator(5);
std::uniform_real_distribution<>	unif(0.0,1.0);
std::normal_distribution<> 		norm{0,0.1};
std::uniform_int_distribution<int> uni_int(0,100000);


double uniformDistribution( double max = 1.0 ){
	std::uniform_real_distribution<>  uniform(0.0,max);
	return(uniform(generator));
}

Vector normalDistributionVector( double std = 0.01 ){
    return(Vector(norm(generator),norm(generator),norm(generator)));
}

// Uniformly distributed point with an option to shrink the window size (for estimation)
Point uniformDistributionPoint(double window = 1.0){
    double pad = (1-window)/2;
    Point p = Point(uniformDistribution(window) + pad,uniformDistribution(window) + pad,uniformDistribution(window) + pad);
	return(p);
}


// Uniformly distributed weighted point with an option to shrink the window size (for estimation)
Weighted_point uniformDistributionWeightedPoint( double window = 1.0, double w_max = 0.01 ){
	Point p = uniformDistributionPoint(window);
    Weight w;
	if (DELAUNAY) {
        w = 0.0001; 
    }
    else {
        w = uniformDistribution(w_max);
    }
	return(Weighted_point(p,w));
}



class Poisson_Delaunay{
private:
    Rt T;
    int n;

public:
    Poisson_Delaunay() {};
    

    void initialize(int n) {
        std::vector<Point> points;
        std::vector<Weighted_point> weighted_points;
        CGAL::Random_points_in_cube_3<Point,Creator> g(0.5);
        std::copy_n( g, n, std::back_inserter(points));

        for (Point& p: points) {
            Weight w = 0.0001;
            weighted_points.push_back(Weighted_point(p,w));
        }

        T.insert(weighted_points.begin(), weighted_points.end());
    }


    int numberOfCells() const {
        return T.number_of_cells();
    }


    int numberOfPoints() const {
        return T.number_of_vertices();
    }
};






class Gibbs_Delaunay{
private:
	double minimum_face_area;
	double maximum_circumradius;
	double theta;
	double intensity;
	
	double energy = 0;
	int forbidden; // number of forbidden cells/edges detected 

	int number_of_active_points = 0;

	double max_weight;

	Rt T;


public:
	Gibbs_Delaunay( double _minimum_face_area = 0.001, double _maximum_circumradius = 0.4, double _theta = 1.0, double _intensity = 500.0, double _max_weight = 0.01):
		minimum_face_area(_minimum_face_area), maximum_circumradius(_maximum_circumradius), theta(_theta), intensity(_intensity), max_weight(_max_weight), forbidden(0) {  }; 

	int numberOfPoints() const {
		return T.number_of_vertices();	
	}
	

	bool isActive(const Cell_handle& c,  std::string type = "centroid", float unit = 1) const {
        bool active = true;
        if (type == "centroid") {
            active = isWithinUnitBox(CGAL::centroid(T.tetrahedron(c)), unit);
        }
        else if (type == "all") {
            for (int i = 0; i < 4; i++) {
                if (!isWithinUnitBox(T.point(c->vertex(i)))) active = false;
            }
        }
        else if (type == "circumsphere") {
            Tetrahedron t = T.tetrahedron(c);
            double r = circumradius(t);
            Point s = CGAL::circumcenter(t);
            double dist_squared = r*r;
            if (s.x() < 0) { dist_squared -= pow(s.x(),2); }
            else if (s.x() > 1) {dist_squared -= pow(s.x() - 1,2);}
            if (s.y() < 0) { dist_squared -= pow(s.y(),2); }
            else if (s.y() > 1) {dist_squared -= pow(s.y() - 1,2);}
            if (s.z() < 0) { dist_squared -= pow(s.z(),2); }
            else if (s.z() > 1) {dist_squared -= pow(s.z() - 1,2);}
            active = (dist_squared >  0);
        }
        else  {
            throw std::invalid_argument("isActive: type has to be either 'centroid', 'all', or 'circumsphere'");
        }
       
		return active;

	}

    
    

	// TODO: Don't have the same function twice
	double updateEnergy(const Rt::Finite_cells_iterator& begin, const Rt::Finite_cells_iterator& end, bool add = true, bool update = true, bool hardcore = true) {
		double energy_update = 0;

        int cellcount = 0;

		for (Rt::Finite_cells_iterator cell = begin; cell != end; ++cell) {
           
            
			if (T.is_infinite(cell)) { continue; }
			Tetrahedron t = T.tetrahedron(cell);
			if (!isActive(cell, "circumsphere")) { continue; } 
            cellcount ++;

            if (hardcore) {
                if (minimumFaceArea(t) < minimum_face_area){
                    std::cout << "Small face: " << minimumFaceArea(t)  << std::endl;
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
	
	double updateEnergy(const std::vector<Cell_handle>& cells, bool add = true, bool update = true, bool hardcore = true) {
		double energy_update = 0;
		for (const Cell_handle& cell: cells){
			if (T.is_infinite(cell)) { continue; }
			Tetrahedron t = T.tetrahedron(cell);
			if (!isActive(cell, "circumsphere")) { continue; } 


            if (hardcore) {
                if (minimumFaceArea(t) < minimum_face_area){
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


	// double expEnergy(double energy){
	// 	return(exp(-energy));
	// }

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
			int counter = 0;
			for (Rt::Finite_vertices_iterator v = T.finite_vertices_begin(); v != T.finite_vertices_end(); ++v){
				if (counter == random_point){
					vertex = v;
					break;
				}
				counter++;
			}
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


	void initialize(bool from_file = true, std::string filename = "files/regular-grid.txt") {
		if (from_file){
			std::ifstream f(filename);
			f.precision(17);
			f >> T;
			std::cout << "Reading from file. Checking if the tesselation is valid: " << T.is_valid() << std::endl;
		} 
		else
		{
			int number_of_points = pow(19,3); 
			std::vector<Point> pointsA;
            std::vector<Point> pointsB;
			std::vector<Weighted_point> weighted_points;
			pointsA.reserve(number_of_points);
            pointsB.reserve(number_of_points);

			CGAL::points_on_cube_grid_3( 1, number_of_points, std::back_inserter(pointsA), Creator() );
			CGAL::points_on_cube_grid_3( 1, number_of_points, std::back_inserter(pointsB), Creator() );
			for (Point& p: pointsA){
				p += Vector(0.5,0.5,0.5);
				Weight w = 0.0001;
				weighted_points.push_back(Weighted_point(p,w));
			}

            // Calculate the distance by which to move the second grid
            // double d = 1.0/(number_of_points - 1); // ( d = (b-a)/(n-1) * (1/2), here b-a=2 and n = number of points)
            // std::cout << "d = " << d << std::endl;
			// for (Point& p: pointsB){
			// 	p += Vector(0.5 + d,0.5 + d,0.5 + d);
			// 	Weight w = 0.0001;
			// 	weighted_points.push_back(Weighted_point(p,w));
			// }

			T.insert(weighted_points.begin(), weighted_points.end());
		}

		number_of_active_points = numberOfActivePoints();

		updateEnergy(T.finite_cells_begin(), T.finite_cells_end());	
		
		std::cout << "Initialization done. Initial tessellation with " << T.number_of_vertices() << " vertices and " << T.number_of_cells() << " cells created. Initial energy: " << ( forbidden ? -1 : energy ) << std::endl;
		std::cout << "Parameters. Theta: " << theta << ", minimum_face_area: " << minimum_face_area << ", maximum_circumradius: " << maximum_circumradius << ", intensity: " << intensity << ", max_weight: " << max_weight << std::endl; 

		if (forbidden){
			std::cout << "Tesselation is not permissible! Number of forbidden elements: " << forbidden << std::endl;
		}

			
	}


	double realEnergy(){
		return(updateEnergy(T.finite_cells_begin(), T.finite_cells_end(), false, false ));
	}


    // This function is now obsolete: hidden points are disabled
 	// void writeHiddenPoints() const {
	// 	std::cout << "Hidden points: " << std::endl;
	// 	int count = 0;
	// 	for (Rt::Finite_cells_iterator c = T.finite_cells_begin(); c != T.finite_cells_end(); ++c){
	// 		std::vector<Weighted_point> v;
	// 		std::copy(c->hidden_points_begin(), c->hidden_points_end(), std::back_inserter(v));
	// 		for (Weighted_point p: v){
	// 			std::cout << p << std::endl;
	// 			++count;
	// 	  	}
	// 	}
	// 	std::cout << "Number of hidden points: " << count;
	// 	

	// }

    // Checks if the point is present in the cells (checks for weight equality, too)
    bool cellHasPoint(const Rt::Cell_handle& c, const Weighted_point& wp){
        bool result = false;
        for (int i = 0; i < 4; i++){
            if (!T.is_infinite(c->vertex(i)) && (T.point(c->vertex(i)) == wp) && (T.point(c->vertex(i)).weight() == wp.weight())){
                    result = true;
                    break;
            }
        }
        return result;
    
    }
    
    // Extracts set of vertices from conflicting cells
    std::set<Rt::Vertex_handle> findConflictingVertices( std::vector<Rt::Cell_handle>& conflicting_cells ){
           std::set<Rt::Vertex_handle> conflicting_vertices;
           if (!conflicting_cells.empty()) {
                   for (Rt::Cell_handle& c: conflicting_cells){
                           for (int i = 0; i < 4; ++i) {
                                   conflicting_vertices.insert(c->vertex(i));
                           }
                   }
            }
           return conflicting_vertices;
    }
    // Gives the opposite face in a cell as a set of points
    std::set<Rt::Vertex_handle> getOppositeFace(const Rt::Cell_handle& c, int i){
        std::set<Rt::Vertex_handle> v;
        v.insert(c->vertex((i+1)%4));
        v.insert(c->vertex((i+2)%4));
        v.insert(c->vertex((i+3)%4));
        return v;
    }
    
    // Extracts set of vertices from boundary faces
    std::set<Rt::Vertex_handle> findBoundaryVertices(std::vector<Rt::Facet>& boundary_faces){
        std::set<Rt::Vertex_handle> boundary_vertices;
        if (!boundary_faces.empty()){
                for (Rt::Facet& f: boundary_faces){
                    std::set<Rt::Vertex_handle> face_vertices = getOppositeFace(f.first, f.second);
                    boundary_vertices.insert(face_vertices.begin(), face_vertices.end());
                }
        }
        return boundary_vertices;
    }

    // Decides if a point can be added, i.e. 1. Is not already present, 2. Won't delete some existing point, 3. Won't be deleted itself
    bool doesNotConflict(const Weighted_point& p){
        Cell_handle containing_cell = T.locate(p);
        if (cellHasPoint(containing_cell, p)) return false;
        
        std::vector<Cell_handle> conflicting_cells;
        std::vector<Rt::Facet> boundary_faces;

        T.find_conflicts(p, containing_cell, std::back_inserter(boundary_faces), std::back_inserter(conflicting_cells));

        if (conflicting_cells.empty()) return false;
        
            
        std::set<Rt::Vertex_handle> conflicting_vertices = findConflictingVertices(conflicting_cells);
        std::set<Rt::Vertex_handle> boundary_vertices = findBoundaryVertices(boundary_faces);
        
        std::set<Rt::Vertex_handle> weaker_vertices;
        std::set_difference(conflicting_vertices.begin(), conflicting_vertices.end(), boundary_vertices.begin(), boundary_vertices.end(), std::inserter(weaker_vertices, weaker_vertices.end()));

        if (!weaker_vertices.empty()) return false;

        return true;
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

    	// Find newly created cells = cells incident to p
		std::vector<Cell_handle> incident_cells;
  		T.incident_cells(vertex, std::back_inserter(incident_cells));
  		updateEnergy(incident_cells);

		return(vertex);
	}


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
	
	// TODO: Fix the bug where energy changes after a move -- is it not already fixed? Check
	Rt::Vertex_handle move(Rt::Vertex_handle vertex_from, Weighted_point point_to){
		remove(vertex_from);
		Rt::Vertex_handle vertex = add(point_to);
		return(vertex);
	}

	void step( std::ofstream& f){
		double a = unif(generator);
		if (a < 1.0/3.0) {
			if (COUT) std::cout << "B,";            // TODO: Improve the output? E.g. not having two streams, having options for streams (more verbose?)
            if (FOUT) f << "B,";
			int prev_total_n = T.number_of_vertices();
			int n = number_of_active_points;
		    	
            Weighted_point x = uniformDistributionWeightedPoint();
			if (doesNotConflict(x)) {

			    if (COUT) std::cout << x << ", ,";
                if (FOUT) f << x << ", ,";

			    double energy_before = energy;
			    Rt::Vertex_handle added_vertex = add(x);
		            	
			    if (COUT) std::cout << energy_before << "," << ( forbidden ? -1 : energy ) << ","; 
                if (FOUT) f << energy_before << "," << ( forbidden ? -1 : energy ) << ","; 
			    
                // See if we accept the proposal
                double b = unif(generator);
                if (COUT) std::cout << b << ",";
                if (FOUT) f << b << ",";
                // double ratio = ( expEnergy(energy)*intensity / ( (n+1) * expEnergy(energy_before)));
				double ratio =  exp(energy_before - energy) * intensity / (n+1);
                if (COUT) std::cout << ratio << ",";
                if (FOUT) f << ratio << ",";
                if ( forbidden || (b > ratio )){ // If not accepted, roll back
                    T.remove(added_vertex);
                    energy = energy_before;
                    forbidden = 0;
                    if (COUT) std::cout << "0,";
                    if (FOUT) f << "0,";
                    if (COUT) std::cout << T.number_of_vertices() << "," << number_of_active_points << ",";
                    if (FOUT) f << T.number_of_vertices() << "," << number_of_active_points << ",";
                }
                else{
                    if (COUT) std::cout << "1" << ",";
                    if (FOUT) f << "1" << ",";
                    if (COUT) std::cout << T.number_of_vertices() << "," << number_of_active_points << ",";
                    if (FOUT) f << T.number_of_vertices() << "," << number_of_active_points << ",";
                }

			    number_of_active_points += (T.number_of_vertices() - prev_total_n);
            }
		}
		else if ( a > 2.0/3.0) {
			if (COUT) std::cout << "D,";	
            if (FOUT) f << "D,";	
			// Choose a random point to remove
			Rt::Vertex_handle vertex = chooseRandomVertexBetter();
			Weighted_point removed_point = T.point(vertex);
			
			if (COUT) std::cout << removed_point << ", ,";
            if (FOUT) f << removed_point << ", ,";

			int prev_total_n = T.number_of_vertices();
			int n = number_of_active_points;
			double energy_before = energy;	
			remove(vertex);
			

			if (COUT) std::cout << energy_before << "," << ( forbidden ? -1 : energy ) << ","; 
            if (FOUT) f << energy_before << "," << ( forbidden ? -1 : energy ) << ","; 
			
			// See if we accept the pproposal
			double b = unif(generator);
			if (COUT) std::cout << b << ",";
            if (FOUT) f << b << ",";
			//double ratio = n*expEnergy(energy) / (expEnergy(energy_before)*intensity); 
			double ratio = n* exp(energy_before - energy) / intensity;
			if (COUT) std::cout << ratio << ",";
            if (FOUT) f << ratio << ",";
			if ( forbidden || ( b > ratio)) {
				T.insert(removed_point);
				energy = energy_before;
				forbidden = 0;
				if (COUT) std::cout << "0,";
                if (FOUT) f << "0,";
				if (COUT) std::cout << T.number_of_vertices() << "," << number_of_active_points << ",";
                if (FOUT) f << T.number_of_vertices() << "," << number_of_active_points << ",";
			}	
			else{
				if (COUT) std::cout << "1"<< ",";
                if (FOUT) f << "1"<< ",";
				if (COUT) std::cout << T.number_of_vertices() << "," << number_of_active_points << ",";
                if (FOUT) f << T.number_of_vertices() << "," << number_of_active_points << ",";

			}

			number_of_active_points += (T.number_of_vertices() - prev_total_n);

		}
		else {
			if (COUT) std::cout << "M,";	
            if (FOUT) f << "M,";	
			// Choose a random point to move
			Rt::Vertex_handle vertex = chooseRandomVertexBetter();
			Weighted_point old_point = T.point(vertex);
			// Generator a random vector by which the point will move
			Vector move_by = normalDistributionVector(); 
			Weighted_point x( bouncePoint(T.point(vertex).point() + move_by)  , T.point(vertex).weight());
			// if (COUT) std::cout << "Proposing to move the point " << std::endl << T.point(vertex) << " by " << std::endl <<  move_by << " to " << std::endl <<  x;
			
            if (doesNotConflict(x)) {
                if (COUT) std::cout << old_point << "," << x << ",";
                if (FOUT) f << old_point << "," << x << ",";
                int prev_total_n = T.number_of_vertices();
                int n = number_of_active_points;

                double energy_before = energy;
                Rt::Vertex_handle new_vertex = move(vertex, x);

                if (COUT) std::cout << energy_before << "," << ( forbidden ? -1 : energy ) << ","; 
                if (FOUT) f << energy_before << "," << ( forbidden ? -1 : energy ) << ","; 

                
                // See if we accept the proposal
                double b = unif(generator);
                if (COUT) std::cout << b << ",";
                if (FOUT) f << b << ",";
                // double ratio = expEnergy(energy) / expEnergy(energy_before);
				double ratio = exp(energy_before - energy);
                if (COUT) std::cout << ratio << ",";
                if (FOUT) f << ratio << ",";
                if ( forbidden || (b > ratio )) {
                    T.remove(new_vertex);
                    T.insert(old_point);
                    energy = energy_before;
                    forbidden = 0;
                    if (COUT) std::cout << "0,";
                    if (FOUT) f << "0,";
                    if (COUT) std::cout << T.number_of_vertices() << "," << number_of_active_points << ",";
                    if (FOUT) f << T.number_of_vertices() << "," << number_of_active_points << ",";
                }	
                else{
                    if (COUT) std::cout << "1" << ",";
                    if (FOUT) f << "1" << ",";
                    if (COUT) std::cout << T.number_of_vertices() << "," << number_of_active_points << ",";
                    if (FOUT) f << T.number_of_vertices() << "," << number_of_active_points << ",";

                }
                number_of_active_points += (T.number_of_vertices() - prev_total_n);
            }
		}
	}


	int numberOfInactivePoints() const {
		int count = 0;	
		for (Rt::Finite_vertices_iterator v = T.finite_vertices_begin(); v != T.finite_vertices_end(); ++v){
			if (!isWithinUnitBox(T.point(v))) { count++; };	
		}
		return count;
	}

	void iterate(int number_of_iterations, std::string filename ){
        std::cout << "Number of iterations: " << number_of_iterations << std::endl;
		clock_t start;
        std::ofstream f(filename);
        if (COUT) std::cout << "step_no,type,pt,pt_mv,energy,energy_after,b,ratio,accept,no_vrt,no_act_vrt,time" << std::endl;
        if (FOUT) f << "step_no,type,pt,pt_mv,energy,energy_after,b,ratio,accept,no_vrt,no_act_vrt,time" << std::endl;
		for (int k = 1; k <= number_of_iterations; ++k){
		 	start = clock();
		 	if (COUT) std::cout << k << "," ;
            if (FOUT) f << k << "," ;
		 	step(f);
		 	if (COUT) std::cout << (double)(clock()-start)/CLOCKS_PER_SEC ;
            if (FOUT) f << (double)(clock()-start)/CLOCKS_PER_SEC;
		 	// std::cout << " " << number_of_active_points << " " <<  numberOfInactivePoints() ;
		 	// std::cout << " Real energy: " << realEnergy(); 
            if (FOUT) f << std::endl;
            if (COUT) std::cout << std::endl;
		}
		
		std::cout << std::endl << "isValid test: " << T.is_valid() << std::endl;
		std::cout << "forbidden :  " << forbidden << std::endl;
	}



	void writeTessellationToFile(std::string filename) const {
		std::ofstream f(filename);
		f.precision(17);
		f << T;

	}

	
    
    bool isRemovable( const Rt::Vertex_handle & v ) {
        bool is_removable = true;
        Weighted_point p = T.point(v);
        remove(v);
        if (forbidden) { is_removable = false; }
        add(p);
        return is_removable;
    }

    // TODO: Bug - this functon and the counter in analyze() give slightly different results
    int numberOfRemovablePoints() {
		int count = 0;	
		for (Rt::Finite_vertices_iterator v = T.finite_vertices_begin(); v != T.finite_vertices_end(); ++v){
			if ( isWithinUnitBox(T.point(v)) & isRemovable(v) ) { count++; };	
		}
		return count;
    }




    double localEnergy( const Rt::Vertex_handle & v ) {
        Weighted_point p = T.point(v);
        double energy_before_removing = energy;
        remove(v);
        double energy_after_removing;
        if (forbidden) { 
            energy_after_removing = std::numeric_limits<double>::infinity();
        }
        else{
            energy_after_removing = energy;
        }
        add(p);

        return(energy_before_removing - energy_after_removing);
    }


    // TODO: Check for existing point?
    double localEnergy( const Weighted_point & p ) {
        if (!doesNotConflict(p)) { return std::numeric_limits<double>::infinity();  }    // Do not allow conflicting points 
        double energy_before_adding = energy;
        Rt::Vertex_handle v = add(p);
        double energy_after_adding;
        if (forbidden) {
            energy_after_adding = std::numeric_limits<double>::infinity();   
        }
        else {
            energy_after_adding = energy;
        }
        remove(v);

        
        return(energy_after_adding - energy_before_adding);
    }


        


    



    double evaluateThetaEquation(double theta_estimate, const std::vector<double> & local_energy_samples, double constant){
        double sum_value = 0;
        for (double h_i: local_energy_samples){
            sum_value += exp(-theta_estimate * h_i) * ( h_i - constant );
            // std::cout << "Sum value: " << sum_value << " after adding local energy " << h_i << std::endl; 
        }
        return sum_value;
    }


    double evaluateThetaKnownZEquation(double theta_estimate, const std::vector<double> & local_energy_samples, double constant){
        double sum_value = 0;
        for (double h_i: local_energy_samples){
            sum_value += intensity*h_i*exp(-theta_estimate* h_i) - constant;
        }
        return sum_value;
    }


    std::tuple<double,double,double,double>  estimate(int samples_arg, float sampling_window) {
        std::cout << "Estimating using the sampling window " << sampling_window << std::endl; 
        int samples_count = 0;
        std::vector<double> samples;

        // Sample local energy for the integrals
        std::cout << "Sampling" << std::endl;
        while (samples_count < samples_arg){
            std::cout << samples_count << " "; 
            Weighted_point p = uniformDistributionWeightedPoint(sampling_window);
            double local_energy = localEnergy(p) / theta; // Divide by theta to obtain energy with theta = 1
            if (std::isfinite(local_energy)) { 
                ++samples_count;
                samples.push_back(local_energy);
            }
        }

        
        /// Theta with z unknown
        // Evaluate the sum over removable points
        // Get their count while iterating over them instead of using numberOfRemovablePoints()
		int number_of_removable_points = 0;	
        std::vector<double> local_energy_removable_points;
		for (Rt::Finite_vertices_iterator v = T.finite_vertices_begin(); v != T.finite_vertices_end(); ++v){
			if (isWithinUnitBox(T.point(v), sampling_window)) { 
                double local_energy = localEnergy(v) / theta; // Divide by theta to obtain energy with theta = 1
                if (std::isfinite(local_energy)) {
                    ++number_of_removable_points;
                    local_energy_removable_points.push_back(local_energy);
                }
            }	
		} 
        double constant_unnormalized = std::accumulate( local_energy_removable_points.begin(), local_energy_removable_points.end(), 0.0);
        std::cout << "Constant unnormalized: " << constant_unnormalized << std::endl;
        double constant = constant_unnormalized / number_of_removable_points;
        std::cout << "Constant normalized: " << constant << std::endl;
     
        
        // std::vector<double> ticks; 
        // double min = -100;
        // double max = 100;
        // int n = 200;
        // double step = (max-min)/n;
        // for (int i = 0; i < n; ++i) {
        //     ticks.push_back(min + i*step);
        // }

        // for (double tick: ticks){
        //     std::cout << evaluateThetaEquation(tick,samples,constant) << std::endl;
        // }
       
        std::cout << std::endl;
        std::cout << "Number of removable points: " << number_of_removable_points << std::endl;


        // The best equation solver
        double error = 1;
        double upper = theta + 15;
        double lower = theta - 15;
        double estimate;
        assert( sign(evaluateThetaEquation(lower,samples,constant)) != sign(evaluateThetaEquation(upper,samples,constant)) ); 

        while (error > 0.00001){
           estimate = (lower + upper) / 2.0;
           error = fabs(evaluateThetaEquation(estimate,samples,constant));
           std::cout << lower << " " <<  upper << " " <<  estimate << " " <<  error << " " << evaluateThetaEquation(lower,samples,constant) << " " << evaluateThetaEquation(upper,samples,constant) << std::endl;
           if (sign(evaluateThetaEquation(lower,samples,constant)) == sign(evaluateThetaEquation(estimate,samples,constant))) {
               lower = estimate;
           }
           else {
               upper = estimate;
           }
        }           


        double theta_estimate = estimate;


        /// Estimate z (intensity)
        // Estimate the integral
        double integral_estimate = 0;
        for (double h_i: samples){
            integral_estimate += exp(-theta_estimate*h_i);
        }
        integral_estimate =  integral_estimate / samples_count;
        
        double z_estimate = number_of_removable_points / integral_estimate;


        
        /// Estimate z with theta known (intensity)
        // Estimate the integral
        integral_estimate = 0;
        for (double h_i: samples){
            integral_estimate += exp(-theta*h_i);
        }
        integral_estimate =  integral_estimate / samples_count;
        double z_known_theta_estimate = number_of_removable_points / integral_estimate;



        /// Theta with z known
        // The best equation solver, again
        error = 1;
        upper = theta + 15;
        lower = theta - 15;
        constant = constant_unnormalized;
        assert( sign(evaluateThetaKnownZEquation(lower,samples,constant)) != sign(evaluateThetaKnownZEquation(upper,samples,constant)) ); 

        while (error > 0.0001){
           estimate = (lower + upper) / 2.0;
           error = fabs(evaluateThetaKnownZEquation(estimate,samples,constant));
           // std::cout << lower << " " <<  upper << " " <<  estimate << " " <<  error << " " << evaluateThetaKnownZEquation(lower,samples,constant) << " " << evaluateThetaKnownZEquation(upper,samples,constant) << std::endl;

           if (sign(evaluateThetaKnownZEquation(lower,samples,constant)) == sign(evaluateThetaKnownZEquation(estimate,samples,constant))) {
               lower = estimate;
           }
           else {
               upper = estimate;
           }
        }           
        double theta_known_z_estimate = estimate;

        return std::make_tuple(theta_estimate,z_estimate,theta_known_z_estimate, z_known_theta_estimate);
    }

    
    // void writeMetadata( std::string filename ) const {
    //     std::ofstream f(filename);
    //     f << "Minimum edge length: " << minimum_edge_length << std::endl;
    //     f << "Maximum circumradius: " << maximum_circumradius << std::endl;
    //     f << "Theta: " << theta << std::endl;
    //     f << "Intensity: " << intensity << std::endl;
    //     f << "Max weight: " << max_weight << std::endl;
    // }
   
   
	void analyze( std::string filename, int samples_arg  )  {
        std::cout << std::endl;
        std::cout << "Analyzing.." << std::endl;

		// Get only active cells
        // Gather data from active cells
		std::vector<Tetrahedron> active_tetrahedra;
		std::set<Face> active_faces;
		std::set<Edge> active_edges;
		std::set<Rt::Vertex_handle> active_vertices;
        std::set<Rt::Vertex_handle> vertices_in_unit_box;

		for (Rt::Finite_cells_iterator cell = T.finite_cells_begin(); cell != T.finite_cells_end(); ++cell) {
            if (isActive(cell,"all")) { 
				Tetrahedron t = T.tetrahedron(cell);				

				active_tetrahedra.push_back(t); 
				for (int i = 0; i < 4 ; ++i) {
					active_vertices.insert(cell->vertex(i));
                    if (isWithinUnitBox(T.point(cell->vertex(i)))) { vertices_in_unit_box.insert(cell->vertex(i)); }
					active_faces.insert(Face(t.vertex(i), t.vertex(i+1), t.vertex(i+2)));
					for (int j = i + 1; j < 4; ++j){
						active_edges.insert(Edge(t.vertex(i),t.vertex(j)));	
					}
				}
				
			} 
		}
		
        // TODO: For some inexplicable reason there are more circumradii than active cells (in Python)
		std::vector<double> tetrahedra_volumes;
        std::vector<double> tetrahedra_surface;
        std::vector<double> tetrahedra_circumradii;
		std::vector<double> face_surfaces;
		std::vector<double> edge_lengths;
		std::vector<double> point_weights;
		std::vector<int> point_degrees;
		

		// Data to output
		// Distributions
        int number_of_removable_points = 0;

		for(Tetrahedron t: active_tetrahedra) { 
            tetrahedra_volumes.push_back(t.volume());
            tetrahedra_circumradii.push_back(circumradius(t));
            tetrahedra_surface.push_back(surfaceArea(t));
        }
		for(Face f: active_faces) { face_surfaces.push_back(faceArea(f)); }
		for(Edge e: active_edges) { edge_lengths.push_back(edgeLength(e)); }
		for(Rt::Vertex_handle v: vertices_in_unit_box) { 
			point_weights.push_back(T.point(v).weight());
			point_degrees.push_back(T.degree(v));
            if (isRemovable(v)) { number_of_removable_points++; } 

		}

		// Scalars
		int number_of_vertices = point_weights.size();
		int number_of_cells = active_tetrahedra.size();

		std::cout << "Number of tetrahedra: " << tetrahedra_volumes.size() << " " << tetrahedra_surface.size() << " " << std::endl;
        std::cout << "No. of faces: " << face_surfaces.size() << std::endl; 
        std::cout << "No. of edges: "  << edge_lengths.size() << std::endl;
        std::cout << "No. of points: " << point_degrees.size() << " " << point_weights.size() << std::endl;


        std::cout << "Point degrees: " << std::endl;
		for(Rt::Vertex_handle v: vertices_in_unit_box) { 
            std::cout << T.degree(v) << " " << T.point(v) << std::endl;
		}




        // Estimate hardcore parameters
        double min_edge_est = *std::min_element( edge_lengths.begin(), edge_lengths.end() );
        double min_face_est = *std::min_element( face_surfaces.begin(), face_surfaces.end() );
        double max_circumradius_est = *std::max_element( tetrahedra_circumradii.begin(), tetrahedra_circumradii.end()  );

        std::cout << std::endl;
        std::cout << "Hardcore estimates" << std::endl;
        std::cout << "Minimum edge: " << min_edge_est << std::endl;
        std::cout << "Minimum face: " << min_face_est << std::endl; 
        std::cout << "Maximum circumradius: " << max_circumradius_est << std::endl;

        

        // Estimate smooth interaction parameters
        std::tuple<double,double,double,double> smooth_estimates  = estimate(samples_arg, 1.0);

        // Output to a file
        std::ofstream f(filename);
        f << "epsilon;" << "alpha;" << "theta;" << "z;" << "max_weight;" << "energy;" << "tetra_volume;" << "tetra_surface;" << "tetra_circum;" <<  "face_surf;" << "edge_length;" << "point_weight;" << "point_degree;" << "cells;" << "vertices;" << "removable;" << "epsilon_est;" << "face_est;" << "alpha_est;" << "theta_est;" << "z_est;" << "theta_known_z_est;" << "z_known_theta_est" << std::endl; 

        f << minimum_face_area << ";" << maximum_circumradius << ";" << theta << ";" << intensity << ";" << max_weight << ";" << energy << ";";
        f << tetrahedra_volumes << ";" << tetrahedra_surface << ";" << tetrahedra_circumradii << ";" << face_surfaces << ";" << edge_lengths << ";" << point_weights << ";" << point_degrees << ";";
        f << number_of_cells << ";" << number_of_vertices << ";" << number_of_removable_points << ";";
        f << min_edge_est << ";" << min_face_est << ";" << max_circumradius_est << ";" << std::get<0>(smooth_estimates) << ";" << std::get<1>(smooth_estimates) << ";" << std::get<2>(smooth_estimates) << ";" << std::get<3>(smooth_estimates);
        f << std::endl;


        std::cout << std::endl;
        std::cout << "Smooth parameter estimates" << std::endl;
        std::cout << "Theta: " << std::get<0>(smooth_estimates) << std::endl;
        std::cout << "z: " << std::get<1>(smooth_estimates) << std::endl;
        std::cout << "Theta (known z): " << std::get<2>(smooth_estimates) << std::endl;
        std::cout << "z (known theta): " << std::get<3>(smooth_estimates) << std::endl;

	}




};
// TODO: Suggesting only addable points
// 

// TODO: Improve arguments, have defaults or so
// Arguments
// 1 Coef
// 2 Exponent
// 3 Theta
// 4 z
// 5 Min face area
// 6 Max circumradius
// 7 Samples count


int main(int agrc, char* argv[]) {
    // Get a timestamp for the files
    // TODO: Improve the filenames / folders
    std::chrono::time_point<std::chrono::system_clock> time_now = std::chrono::system_clock::now();
    std::time_t time_now_t = std::chrono::system_clock::to_time_t(time_now);
    std::tm now_tm = *std::localtime(&time_now_t);
    char buf[512];
    std::strftime(buf, 512, "_%Y%m%d_%H_%M_%S", &now_tm);
    
    // std::chrono::milliseconds ms = std::chrono::duration_cast< std::chrono::milliseconds >(std::chrono::system_clock::now().time_since_epoch());
    // std::string ms_s = std::to_string(ms.count()).substr(10,12);


    int random_int = uni_int(generator);

    int coef = std::stoi(argv[1]);
    int expon = std::stoi(argv[2]);
    std::string filename(buf);
    filename = "_" + std::to_string(coef) + "_" + std::to_string(expon) + filename + "_" + std::to_string(random_int);

    double theta = std::stod(argv[3]);
    double z = std::stod(argv[4]);
    double min_face = std::stod(argv[5]);
    double max_circum = std::stod(argv[6]);
    // Min face, max circum, theta
    Gibbs_Delaunay GD(min_face, max_circum, theta, z);

    if (ANALYZE) {
        GD.initialize(true, "files/gibbs.txt");

        std::cout << "Number of points, total: " << GD.numberOfPoints() << std::endl;
        std::cout << "Number of active points (within unit box):  " << GD.numberOfActivePoints() << std::endl;

        GD.analyze( "files/cell_data" + filename + ".txt" , std::stoi(argv[8]));
    }
    else {
        // GD.initialize(true, "files/gibbs.txt");
        // GD.initialize(true, "files/regular-grid.txt");  
        GD.initialize(false);
        GD.iterate(coef*pow(10,expon), "files/log" + filename + ".csv");

        std::cout << "Number of points, total: " << GD.numberOfPoints() << std::endl;
        std::cout << "Number of active points (within unit box):  " << GD.numberOfActivePoints() << std::endl;

        GD.writeTessellationToFile("files/gibbs" + filename + ".txt");
        GD.analyze( "files/cell_data" + filename + ".txt" , std::stoi(argv[8]));
    }


    // Running Poisson
    // for (int i =1; i <= 100; i++){
    //     for (int j = 1; j <= 100; j++){
    //         Poisson_Delaunay PD;
    //         PD.initialize(10*i);
    //         std::cout<< PD.numberOfPoints() << " " << PD.numberOfCells() << std::endl;
    //     }
    // }

}
