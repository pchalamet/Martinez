#include "polygon.h"
#include "utilities.h"
#include "martinez.h"
#include "connector.h"
#include "greiner.h"
#include "gpc.h"
#include <fstream>

int main (int argc, char* argv[])
{
	if (argc < 4) {
		std::cerr << "Syntax: " << argv[0] << " subject_pol clipping_pol result_pol [I|U|D|X]\n";
		return 1;
	}
	if (argc > 4 && argv[4][0] != 'I' && argv[4][0] != 'U' && argv[4][0] != 'D' && argv[4][0] != 'X') {
		std::cerr << "Syntax: " << argv[0] << " subject_pol clipping_pol result_pol [I|U|D|X]\n";
		std::cerr << "The last parameter is optional. It is a character. It can be I (Intersection), U (Union), D (Difference) or X (eXclusive or)\n";
		return 2;
	}
	Martinez::BoolOpType op = Martinez::INTERSECTION;
	gpc_op opVatti = GPC_INT;
	if (argc > 4) {
		switch (argv[4][0]) {
			case 'I':
				op = Martinez::INTERSECTION;
				opVatti = GPC_INT;
				break;
			case 'U':
				op = Martinez::UNION;
				opVatti = GPC_UNION;
				break;
			case 'D':
				op = Martinez::DIFFERENCE;
				opVatti = GPC_DIFF;
				break;
			case 'X':
				op = Martinez::XOR;
				opVatti = GPC_XOR;
				break;
		}
	}

	int ntests = 0; // number of tests
	Polygon subj (argv[1]);
	Polygon clip (argv[2]);
	Polygon martinezResult;
	clock_t start;
	int GreinerResult;
	float Martacum = 0;
	float Greineracum = 0;
	float Vattiacum = 0;
	while (Martacum < 1.0f && Greineracum < 1.0f && Vattiacum < 1.0f) {
		ntests++;
		martinezResult.clear ();
		// Martínez-Rueda's algorithm
		start = clock ();
		Martinez mr (subj, clip);
		mr.compute (op, martinezResult);
		Martacum += clock () - start;
		// Greiner-Hormann's algorithm
		Polygon greinerResult;
		start = clock ();
		GreinerHormann gh (subj, clip);
		GreinerResult = gh.boolop (op, greinerResult);
		Greineracum += clock () - start;
		// Vatti's algorithm
		gpc_polygon subject, clipping, result;
		gpc_set_polygon (subj, &subject);
		gpc_set_polygon (clip, &clipping);
		start = clock ();
		gpc_polygon_clip (opVatti, &subject, &clipping, &result);
		Vattiacum += clock () - start;
		gpc_free_polygon (&subject);
		gpc_free_polygon (&clipping);
		gpc_free_polygon (&result);
	}
	std::cout << "Martínez-Rueda's time: " << Martacum / double (CLOCKS_PER_SEC) / ntests << std::endl;
	std::cout << "Vatti's time: " << Vattiacum / double (CLOCKS_PER_SEC) / ntests << std::endl;
	switch (GreinerResult) {
		case -1:
			std::cout << "Sorry, the Greiner-Hormann's method needs perturbation, and it is not implemented." << std::endl;
			break;
		case -2:
			std::cout << "Sorry, the Greiner-Hormann's method cannot work with this operation and polygons with more than one region." << std::endl;
			break;
		default:
			std::cout << "Greiner-Hormann's time: " << Greineracum / double (CLOCKS_PER_SEC) / ntests << std::endl;
			break;
	}
	std::cout << "Possible intersections: " << subj.nvertices () << " x " << clip.nvertices () << " = " << subj.nvertices()*clip.nvertices() << std::endl;
	std::cout << "Number of tests: " << ntests << std::endl;
	ofstream f (argv[3]);
	if (!f) 
		std::cerr << "can't open " << argv[3] << '\n';
	else
		f << martinezResult;
	return 0;
}

