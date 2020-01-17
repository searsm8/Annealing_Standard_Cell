//anneal_layout.cpp
//Mark Sears
//implement layout algorithm with simulated annealing
//layout determined by slicing method

#include <stdio.h>
#include <stdlib.h>  
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <string>
#include <list>
#include <map> 
#include <vector> 
#include <iterator>
#include <iostream>
#include <ctime>
#include <stack>
#include <numeric>

using namespace std;

class BLOCK;
class NET;

//***************************//
//  DECLARE GLOBAL VARIABLES //
//***************************//

int move_weights[] = {0, 70, 20, 10};

float REDUCTION_FACTOR = .97;  //factor that reduces annealing temperature
int     STEPS_PER_TEMP = 3000; //number of annealing steps performed at each temperature

float Si = 4; 	//maximum allowed Aspect Ratio height / width
float Ri = .25; //minimum allowed Aspect Ratio
int num_shapes = 5; //number of shapes a block can form using different A.R.'s

vector <BLOCK*> blocks; //master set of all cells in the design
vector <string> ops; //strings that contains H and V cut operations
vector <NET*> nets; //all nets connecting blocks in the design

BLOCK* best_layout;
vector <BLOCK*> best_blocks; //best block vector found so far
vector <string> best_ops;    //best ops vector found so far
float best_cost;

vector <int> previous_ops_count; //keeps a cumulative count of operators. Used to ensure balloting property


int num_active_blocks = 0; //count how many blocks are in memory....used for debugging
int num_total_blocks = 0; //count how many blocks have every been made.

//used to track program execution time
std::clock_t start_time = std::clock();
std::clock_t timestamp = std::clock();
std::clock_t time_ref = std::clock();
double bestArea_time = 0; //track how long is spent in Stockmeyer
double solidifyBlocks_time = 0; //track how long is spent in Stockmeyer

//a "BLOCK" represents a cell or combination of cells
//each block may contain smaller blocks, which may contain smaller blocks...
//this forms a "tree" of blocks, with the top block representing the whole layout
class BLOCK
{
	public:
	
	string name;
	
	float x; //x and y coordinates of the bottom left corner
	float y;
	float width;
	float height;
	float area;			
	
	bool cut_type; // 0 for H cut, 1 for V cut
	float max_AR; //Aspect Ratio limits
	float min_AR;	
	
	vector < pair<float, float> > shapes; //pair<width, height>
						//set of possible shapes that are within A.R. limits
						//should be sorted by ascending widths
						
	vector <BLOCK*> sub_blocks; //set of sub-blocks that make up this large block (usually 2)
				//base cells have an empty set of sub-blocks
	
	//constructor
	
	BLOCK()
	{	
		num_active_blocks++;
		num_total_blocks++;
		name = "block_" + to_string((long long int)num_total_blocks);
		x = 0;
		y = 0;
		width = 0;
		height = 0;
		area = 0;
		cut_type = false;
		max_AR = Si;
		min_AR = Ri;
		
		
		//printf("Constructor for new BLOCK: %s\n", name.c_str());
	}
	
	/*			
	//destructor	
	~BLOCK()
	{
		num_active_blocks--;
	
	//printf("Deconstructor called for BLOCK: %s\tnum_total_blocks: %d\tnum_active_blocks: %d\t", name.c_str(), num_total_blocks, num_active_blocks);
	
//USE free() HERE??S?S?S?
		//printf("b addr: %p\t", this);
		
		//printf("x addr: %p\t", &x);
		//delete &x; x = NULL; 
		//delete &y; y = NULL;
		//delete &shapes; shapes = NULL;
		
		//printf("\n");
	}
	*/
						
};

class NET {
	public:
	string name;
	
	float wire_length; //estimated length of wire needed for this net
	
	vector <BLOCK*> blocks; //pointers to all blocks on this net			
};

//***************************//
//     FUNCTION HEADERS	     //
//***************************//

void readBlocksFile(string benchmark_path);
float estimateTotalWirelength(BLOCK* b);
float estimateWirelength(NET* net);
float printNet(NET* net);
BLOCK* combineBlocksVCut(BLOCK* b1, BLOCK* b2);
BLOCK* combineBlocksHCut(BLOCK* b1, BLOCK* b2);
void printBlock(BLOCK* b);
BLOCK* findBestArea();
void solidifyBlocks(BLOCK* b);
int pickMove();
void performMove(float temp);
float combinedCost(BLOCK* b);
void move1();
void move2();
void move3();
void update_previous_ops_count();
void initializeOps();
void initializeAnnealing();
float initializeTemp();
void performAnnealingStep();
bool evaluateMove(float temp);
void acceptMove(BLOCK* new_layout, float new_cost);
void rejectMove(BLOCK* new_layout);
void deallocateBlock(BLOCK* b);
void printAnnealingResults(int trial_number, string benchmark_name, float initial_wirelength);
void printTimestamp();
void printBlockShapes(BLOCK* b);
float bestPossibleArea();
void exportBlockShapes(BLOCK* b, char* output_file_name);
void printBlockToFile(BLOCK* b, FILE* output_file);
void anneal_slicing(string benchmark_path);



void readBlocksFile(string benchmark_path)
{	

	//reset global data structures
	blocks.clear();
	nets.clear();
	
	FILE *blocks_file; 
	blocks_file = fopen(benchmark_path.c_str(), "r");	
	if (blocks_file == NULL) { fprintf(stderr, "Can't open blocks_file: %s!\n", benchmark_path.c_str()); exit(1); }
	
	printf("\nReading fpga_benchmarks file: %s\n", benchmark_path.c_str());
			
	char new_block_area [100];
	char new_block_min_width [100];
	char new_block_min_height [100];
	char fpga_name [100];
	int num_blocks;
	int num_wires;
					
	//fscanf(blocks_file, "%s %s %s %s\n", &new_block_area, &new_block_min_AR, &new_block_max_AR);
	
	fscanf(blocks_file, "%s\n", &fpga_name);//read in fpga name
	printf("FGPA name: %s\n", fpga_name);
	
	fscanf(blocks_file, "%d\n", &num_blocks);//read in number of blocks
	printf("Reading blocks...\n");
	
	//read in all blocks
	for(int i = 0; i < num_blocks; i++)
	{		
				
		BLOCK* new_block = new BLOCK;
		//BLOCK* new_block;
		
		fscanf(blocks_file, "%s %s %s\n", &new_block_area, &new_block_min_width, &new_block_min_height);
		
		//convert strings to floats		
		new_block->area = stof(new_block_area);
		new_block->name = "b" + to_string((long long int)i);
		float min_width = stof(new_block_min_width);	
		float min_height = stof(new_block_min_height);
		if(min_width  <= 0) min_width = 1;
		if(min_height <= 0) min_height = 1;
		
		
		
		//limited A.R. (0.25 to 4)
		new_block->min_AR = Ri; new_block->max_AR = Si;
		
		
		
		
		//relaxed A.R. (very low to very high)
		//new_block->min_AR = (min_height*min_height) / new_block->area; new_block->max_AR = new_block->area / (min_width*min_width);
		
		
		
//printf("%s: Area = %.1f \t Min AR = %.4f \t Max AR = %.4f\n", new_block->name.c_str(), new_block->area, new_block->min_AR, new_block->max_AR);
							
		int divisor = (num_shapes-1)/2;
//printf("%s: \n", new_block->name.c_str());
		//generate shapes with different AR for the given area
		int shapes_used = 0;
		float last_w;
		
		//A.R. < 1
		for(float i = new_block->min_AR; i < 1; i += (1-new_block->min_AR)/divisor)
		{
		//printf("%d ", ++shapes_used);
			pair<float, float> new_shape;
			new_shape.first = sqrt(new_block->area * i) ;
			new_shape.second = new_block->area / new_shape.first;
			new_block->shapes.push_back(new_shape);
			last_w = new_shape.first;
		//printf("<%f, %f>\n", new_shape.first, new_shape.second);
		}
		
		//A.R. >= 1
		for(float i = 1; i <= new_block->max_AR; i += (new_block->max_AR-1)/divisor)
		{
		//printf("%d ", ++shapes_used);
			pair<float, float> new_shape;
			new_shape.first = sqrt(new_block->area * i) ;
			new_shape.second = new_block->area / new_shape.first;
			if(new_shape.first - last_w > .001)
			{
			
			new_block->shapes.push_back(new_shape);
		//printf("#<%f, %f> %f\n", new_shape.first, new_shape.second, new_shape.first - last_w);
			}
		}
	//printf("\n");	
		blocks.push_back(new_block);
	}
	
	printf("Reading wires...\n");
	
	fscanf(blocks_file, "\n%d", &num_wires);//read in number of wires
	//printf("\nnum_wires: %d\n", num_wires);
	
	//read in all wires
	for(int i = 0; i < num_wires; i++)
	{		
		int num_blocks_on_net;
		fscanf(blocks_file, "\n%d\n", &num_blocks_on_net);		
		
		NET* new_net = new NET;
		//NET* new_net;
		
		new_net->name = "n" + to_string((long long int)i);
		
				
		for(int j = 0; j < num_blocks_on_net; j++)
		{
			int block_index;
			fscanf(blocks_file, "%d ", &block_index);		
			new_net->blocks.push_back(blocks[block_index]);
		}
				
		nets.push_back(new_net);			
	}
	
	//printf("\ndone!\n");
	return;
}

double Variance(std::vector<double> samples)
{
     int size = samples.size();

     double variance = 0;
     double t = samples[0];
     for (int i = 1; i < size; i++)
     {
          t += samples[i];
          double diff = ((i + 1) * samples[i]) - t;
          variance += (diff * diff) / ((i + 1.0) *i);
     }

     return variance / (size - 1);
}

double StandardDeviation(std::vector<double> samples)
{
     return sqrt(Variance(samples));
}

float estimateTotalWirelength(BLOCK* b)
{
	b->name = "layout";
	solidifyBlocks(b); //must set x and y coordinates of each block before estimating wirelen!

	float total_wirelength = 0;
	
	for(int i = 0; i < nets.size();i++)
	{
		total_wirelength += estimateWirelength(nets[i]);		
	}
//printf("\ntotal_wirelength = %.1f", total_wirelength);
	return total_wirelength;
}

//estimate the wire needed to connect a net using bounding box (half perimeter) method
float estimateWirelength(NET* net)
{
//find the min/max x and y for this net
	float min_x = net->blocks[0]->x;
	float min_y = net->blocks[0]->y;
	float max_x = net->blocks[0]->x;
	float max_y = net->blocks[0]->y;
		
	for(int i = 1; i < net->blocks.size(); i++)
	{
		if(net->blocks[i]->x < min_x)
			min_x = net->blocks[i]->x;
		if(net->blocks[i]->y < min_y)
			min_y = net->blocks[i]->y;
			
		if(net->blocks[i]->x > max_x)
			max_x = net->blocks[i]->x;
		if(net->blocks[i]->y > max_y)
			max_y = net->blocks[i]->y;
	}
	
	return (max_x - min_x) + (max_y - min_y);
}

float printNet(NET* net)
{
	printf("\nNET %s:", net->name.c_str());
	
	for(int i = 0; i < net->blocks.size(); i++)
	{
		printf(" %s", net->blocks[i]->name.c_str());
	}		
}

int shapes_carried_forward = 0;
vector <int> all_shapes_carried_forward;
bool print_carried_num = false;

//combineBlocksVCut (Vertical Cut) takes in two blocks with known areas and Aspect Ratio limits
//and finds the potential best combinations of those blocks for potentially minimum area
//returns a pointer to a new block that has the combined generated shapes
BLOCK* combineBlocksVCut(BLOCK* b1, BLOCK* b2)
{
	BLOCK* combined_block = new BLOCK;
	//BLOCK* combined_block;
	
	int index_1 = 0; //index for block1
	int index_2 = 0; //index for block2
	
	//read the possible shapes in ascending widths
	
	while(index_1 < b1->shapes.size() && index_2 < b2->shapes.size())
	{
		//width of combined block is the sum of widths
		float new_width = b1->shapes[index_1].first + b2->shapes[index_2].first;
		float new_height;
		
		//height of combined block is the larger height
		if(b1->shapes[index_1].second > b2->shapes[index_2].second) // if height of b1 is larger
		{
			new_height = b1->shapes[index_1].second; 
			index_1++;
		}
		else // if height of b2 is larger
		{
			new_height = b2->shapes[index_2].second;
			index_2++;
		}
		
		combined_block->shapes.push_back(make_pair(new_width, new_height));		
		//these shapes must be sorted in ASCENDING widths
	}
shapes_carried_forward += combined_block->shapes.size();

	combined_block->sub_blocks.push_back(b1); //initialize sub_blocks
	combined_block->sub_blocks.push_back(b2);
	combined_block->cut_type = 1; // signifies V cut
	//combined_block->name = b1->name + " V " + b2->name;	
	return combined_block;
}

//combineBlocksHCut (Horizontal Cut)
BLOCK* combineBlocksHCut(BLOCK* b1, BLOCK* b2)
{
	BLOCK* combined_block = new BLOCK; //assign blocks to the heap
	//BLOCK* combined_block;
	
	int index_1 = b1->shapes.size() - 1; //index for block1
	int index_2 = b2->shapes.size() - 1; //index for block2
	
	//read the shapes in ascending heights	
	while(index_1 >= 0 && index_2 >= 0)
	{
		//height of combined block is the sum of heights
		float new_height = b1->shapes[index_1].second + b2->shapes[index_2].second;
		float new_width;
		
		//width of combined block is the larger width
		if(b1->shapes[index_1].first > b2->shapes[index_2].first) // if width of b1 is larger
		{
			new_width = b1->shapes[index_1].first; 
			index_1--;
		}
		else // if height of b2 is larger
		{
			new_width = b2->shapes[index_2].first;
			index_2--;
		}
		
		//combined_block->shapes.push_back(make_pair(new_width, new_height));
		combined_block->shapes.insert(combined_block->shapes.begin(), make_pair(new_width, new_height));
		//these shapes must be sorted in ASCENDING widths
	}
if(print_carried_num) printf("H# carried forward: %d\n", combined_block->shapes.size());
shapes_carried_forward += combined_block->shapes.size();

	combined_block->sub_blocks.push_back(b1); //initialize sub_blocks
	combined_block->sub_blocks.push_back(b2);	
	combined_block->cut_type = 0; // signifies H cut
	//combined_block->name = b1->name + " H " + b2->name;
	return combined_block;
}


//blockLayout combines blocks recursively until a full layout is complete for a given layout string
BLOCK* findBestArea()
{
	time_ref = std::clock(); //track time spent in this function
	 
	 stack <BLOCK*> block_stack; //use a stack to read the layout string. 
	 shapes_carried_forward = 0;
	 
	 //Combine all blocks based on the operator chains
	 for(int i = 0; i < blocks.size(); i++)
	 {
	 //Put the next cell in a block on the stack		
		block_stack.push(blocks[i]);
		
	 //resolve all operators
	 	for(int op_index = 0; op_index < ops[i].size(); op_index++)
		{
	
	//pop the 2 blocks on top and combine them				
	 	BLOCK* b1 = block_stack.top();
		block_stack.pop();
		BLOCK* b2 = block_stack.top();
		block_stack.pop();
		
	//combine the two blocks based on the operator
	 //then put the new combined block back on the stack
		if(ops[i][op_index] == 'H') //read an op from the chain
			block_stack.push(combineBlocksHCut(b1, b2));
		else 	block_stack.push(combineBlocksVCut(b1, b2));					 

	 	}
	 
	 }
	 		 	 
	 BLOCK* layout = block_stack.top();
	 layout->name = "layout";
	 
	 //from all possible shapes of the full layout, pick the smallest area shape
	 //return this as the area of the full layout
	 for(int i = 0; i < layout->shapes.size(); i++)
	 {
	 	float area = layout->shapes[i].first * layout->shapes[i].second;
		
		if(i == 0 || area < layout->area) //modify area if a smaller area is found
		{
			layout->area = area;
			layout->width = layout->shapes[i].first;
			layout->height = layout->shapes[i].second;	
		}
	 }
	 
	 layout->x = 0;
 	 layout->y = 0; //the full layout block always starts at 0, 0
if(print_carried_num) printf("Total shapes carried forward: %d\n", shapes_carried_forward); 
if(print_carried_num) printf("*******END********\n\n");
all_shapes_carried_forward.push_back(shapes_carried_forward);

bestArea_time += clock() - time_ref;
	 return layout;
}


//Use Stockmeyer Algorithm to "work backwards" through the layout of blocks
//Write precise xy block coordinates for each basic block
//this "solidifies" the blocks into definite shapes and coordinates.
//only solidifies the two blocks that make up the current block,
//then recursively calls solidify again until no lower blocks exist
//Once blocks are solidified, wirelength can be estimated
void solidifyBlocks(BLOCK* b)
{	
	if(b->sub_blocks.size() == 0) 
	{
		return; //if no lower blocks, done!		
	}
	
	//solidify the current block's area
	b->area = b->width * b->height;	

	//solidify w and h for sub blocks
	if(b->cut_type) //if V cut, solidify sub-blocks with same height (or next lesser height)
	{
		for(int i = 0; i < b->sub_blocks.size(); i++)
			for(int j = 0; j < b->sub_blocks[i]->shapes.size(); j++)			
			{
				if(b->height >= b->sub_blocks[i]->shapes[j].second)
				{	
					b->sub_blocks[i]->width  = b->sub_blocks[i]->shapes[j].first;
					b->sub_blocks[i]->height = b->sub_blocks[i]->shapes[j].second;
					break;
				}			
			}
			
	}
	else //if H cut, solidify sub-blocks with same width
	{
		for(int i = 0; i < b->sub_blocks.size(); i++)
			for(int j = b->sub_blocks[i]->shapes.size()-1; j >= 0; j--)			
			{
				if(b->width >= b->sub_blocks[i]->shapes[j].first)
				{	
					b->sub_blocks[i]->width  = b->sub_blocks[i]->shapes[j].first;
					b->sub_blocks[i]->height = b->sub_blocks[i]->shapes[j].second;
					break;
				}				
			}			
	}
	
	//solidify x and y coordinates for sub-blocks
	//first sub-block is always same x,y as the top block
	b->sub_blocks[0]->x = b->x;
	b->sub_blocks[0]->y = b->y;
	
	//if V-cut, set 2nd sub-block to the right
	if(b->cut_type)
	{
		b->sub_blocks[1]->x = b->x + b->sub_blocks[0]->width;
		b->sub_blocks[1]->y = b->y;
	}	
	else //else H-cut, set 2nd sub-block above
	{
		b->sub_blocks[1]->x = b->x;
		b->sub_blocks[1]->y = b->y + b->sub_blocks[0]->height;
	}
	
	//recursively soldify all lower-level sub-blocks
	for(int i = 0; i < b->sub_blocks.size(); i++)
		solidifyBlocks(b->sub_blocks[i]);
			
}//end solidifyBlocks

//returns a number corresponding to a move picked randomly based on weights
//e.g. weight1 = 50, weight2 = 25, weight3 = 25
//will pick move1 at 50%, move2 at 25%
int pickMove()
{
	int sum_of_weights = 0;
	for(int i = 1; i <= 3; i++)
		sum_of_weights += move_weights[i];
		
	int pick = rand() % sum_of_weights;
	
	for(int i = 1; i <= 3; i++)
	{
		if(pick < move_weights[i])
			return i;
		pick -= move_weights[i];
	}
}

void performMove(float temp)
{	
	switch(pickMove())
	{
		case 1: move1(); break;
		case 2: move2(); break;
		case 3: move3(); break;
		default:move1();
	}
	
}

//move M1
//swaps two random adjacent blocks
void move1()
{	

	int i = rand() % (blocks.size()-1);

	BLOCK* temp = blocks[i];
	blocks[i] = blocks[i+1];
	blocks[i+1] = temp;

	return;	
}

//move M2
//complements a random operator chain
void move2()
{
	//choose a random op string until a non-empty one is found (never pick ops[0])
	int i = ( rand() % (ops.size()-1) ) + 1;
	while(ops[i].size() == 0)
		i = ( rand() % (ops.size()-1) ) + 1;

	//complement the H's and V's
	string complement = "";
	for(int j = 0; j < ops[i].size(); j++)
	{
		if(ops[i][j] == 'H')
			complement += 'V';
		else complement += 'H';
	}
	
	ops[i] = complement;

	return;
}


//move M3
//swaps a random cell with an adjacent operator
//must maintain normalized Polish expression and balloting property
//pick a random op. Try to move it forward or backward.
//If both are illegal moves, pick another random op 
void move3()
{
	update_previous_ops_count(); //ensure that the op counts are correct!
	
	//attempt to swap a cell and operator. If move is not legal, pick another random op string
	int attempts = 0; 		
	while(attempts < ops.size()) //try this many times, then give up and cancel move
	{
		//choose a random op string until a non-empty one is found (never pick ops[0])
		int i = ( rand() % (ops.size()-1) ) + 1;
		while(ops[i].size() == 0)
			i = ( rand() % (ops.size()-1) ) + 1;				
		
		
		if(previous_ops_count[i-1] < i-1) //check if this move would preserve balloting
		    { 
			if(ops[i-1].size() == 0 || ops[i][0] != ops[i-1][ops[i-1].size()-1]) //check if this move would preserve normalization
			{	//printf("moving backward!\n");
				ops[i-1].push_back(ops[i][0]);
				ops[i].erase(ops[i].begin());								
				previous_ops_count[i-1]++;
				
				break;			
			}
		    } 
		
		
		if(i != ops.size() - 1) //can't move the last op forward!
		{	//printf("moving forward!\n");		
			if(ops[i][ops[i].size()-1] != ops[i+1][0]) //check if this move would preserve normalization
			{
				ops[i+1].insert(0, 1, ops[i][ops[i].size()-1]);
				ops[i].erase(--ops[i].end());								
				previous_ops_count[i]--;
				
				break;		
			} 
		}
		
		attempts++;		
	}	
	
	return;
} //end move3

//moves through the ops vector, counting how many ops at each index
//useful for preserving the balloting property such that there are never more ops than blocks
void update_previous_ops_count()
{
	for(int i = 1; i < ops.size(); i++)	
	{
		previous_ops_count[i] = previous_ops_count[i-1] + ops[i].size();
	}
}


//initializeOps sets all ops to "H"
//blocks: P0, P1, P2, P3, P4...
//ops:     C0, C1, C2, C3, C4...
void initializeOps()
{
	printf("Initializing ops...\n");
	ops.clear();
	
	ops.push_back(""); //first op must be empty
	previous_ops_count.push_back(0); //initialize	
	
	//initialize all ops to "H"
	for(int i = 1; i < blocks.size(); i++)
	{
		ops.push_back("H");
		previous_ops_count.push_back(previous_ops_count[i-1] + 1);
	}
			
	return;
}

//print information about the blocks and operations to console
void printOps()
{
		printf("\n");
		printf("Area = %.3f \t Width = %.3f \t Height %.3f\n", best_layout->area, best_layout->width, best_layout->height);
		
		printf("best_layout:");
		for(int i = 0; i < best_blocks.size(); i++)
			printf("'%s(%.1f x %.1f)' '%s' ", best_blocks[i]->name.c_str(), best_blocks[i]->width, best_blocks[i]->height), best_ops[i].c_str();
		
		printf("\n********************\n");
		
}


float alpha; //relative cost weight of area
float beta; //relative cost weight of wirelength
//perform steps required before annealing can begin
void initializeAnnealing()
{
	printf("\n\n********************\nBegin Annealing...\n********************\n");
	printf("num_shapes: %d\n", num_shapes);
	printf("REDUCTION_FACTOR: %.3f\n", REDUCTION_FACTOR);
	
	best_blocks = blocks;
	best_ops = ops;
	
	alpha = 0.8;
	beta  = 0.2;
	
	best_layout = findBestArea();
	solidifyBlocks(best_layout);
	best_cost = combinedCost(best_layout);
}

//perform many random moves, record deltas
//set starting temperature equal to 30*std_dev
float initializeTemp()
{
	int num_initialize_moves = 10000;
	
	printf("Initializing temp...(%d random moves)\n", num_initialize_moves);
	vector <double> deltas;
	float new_cost, delta;
		
	//perform random cell swaps
	for(int i = 0; i < num_initialize_moves; i++)
		performMove(100000); //random move at very high temp
		
	//perform random cell swaps and record deltas	
	for(int i = 0; i < num_initialize_moves; i++)
	{
		performMove(100000); //random move at very high temp
		
		BLOCK* new_layout = findBestArea();
				
		new_cost = combinedCost(new_layout);
	
		delta = new_cost - best_cost;
		deltas.push_back(delta);
		
		acceptMove(new_layout, new_cost); //always accept the move during initialization	
	}
	
	double std_dev = StandardDeviation(deltas);
		
	float start_temp = 30*std_dev; //for a temperature of 30*std_dev,
			//there will be a 90% chance to accept a very bad change of 3*std_dev
	printf("Starting Temp: %.3f\n", start_temp);
	
	//initialize alpha and beta based on the current area and wirelength
	//initializeWeights(best_layout);
		
	return start_temp;
}


float combinedCost(BLOCK* b)
{
	return alpha*b->area + beta*estimateTotalWirelength(b);
}


//statistics variables. Allows for analysis of annealing steps
int num_moves_accepted = 0;
int num_moves_accepted_Pr = 0;
int num_moves_rejected = 0;
float accepted_deltas = 0;
float accepted_Pr_deltas = 0;
float rejected_deltas = 0;

bool performAnnealingStep(float temp)
{
	//perform a random move
	performMove(temp);
	
	//evaluate the move
	return evaluateMove(temp);	
}

//evaluateMove compares the new layout's area and wirelength to the previous.
//then accepts or rejects the move 
bool evaluateMove(float temp)
{
	BLOCK* new_layout = findBestArea();	
	
	float new_cost = combinedCost(new_layout);
	//printf("Estimated Wirelength: %.3f\n", new_wirelength);
	
	//if a better layout was found, keep it!
	//float delta = new_layout->area - best_layout->area;
	float delta = new_cost - best_cost;
	 		
	if(delta <= 0) //if an improvement has been made
	{		
		acceptMove(new_layout, new_cost);
		
		num_moves_accepted++;
		accepted_deltas += delta;					
		return true; //a move was made	
	}

	else //else no imrovement was made. keep the new layout conditionally
	{
		double Pr_rejecting = exp((float)(-1*delta)/(float)temp);
		
		int p = 100000;		
		if((float)(rand()%p) / (float)p < Pr_rejecting) //chance to accept or reject
		{
			acceptMove(new_layout, new_cost);

			num_moves_accepted_Pr++;
			accepted_Pr_deltas += delta;	
			return true; //a move was made
		}
		else //else reject the bad move, revert to best layout
		{			
			rejectMove(new_layout);
		
			num_moves_rejected++;
			rejected_deltas += delta;			
			return false; //no move was made
		}
		
	}	
}//end evaluateMove()


void acceptMove(BLOCK* new_layout, float new_cost)
{
//printf("acceptMove() for block %s\n", new_layout->name.c_str());
	//deallocate memory of the old best layout
	deallocateBlock(best_layout);
	
	//set new best layout
	best_layout = new_layout;		
	best_blocks = blocks;
	best_ops = ops;
	best_cost = new_cost;	
		
}

void rejectMove(BLOCK* new_layout)
{
//printf("rejectMove()\n");
	//deallocate memory for the rejected layout
	deallocateBlock(new_layout);
	
	blocks = best_blocks;
	ops = best_ops;	
}

void deallocateBlock(BLOCK* b)
{
	//return;
	
	if(b->sub_blocks.size() == 0) return; //don't delete the base level blocks!!!
	
	//recursively deallocate all sub-blocks
	for(int i = 0; i < b->sub_blocks.size(); i++)
		deallocateBlock(b->sub_blocks[i]);
	
	//deallocate the memory for this block
	delete b; b = NULL;
			
	num_active_blocks--;
	
	/*
	for(int i = 0; i < b->sub_blocks.size(); i++)
	{
		free(b->sub_blocks[i]);
		b->sub_blocks[i] = NULL; //assign to NULL to avoid double free
	}
	*/
	
	//free(b);
	
	return;
}

void resetTimestamp()
{
	timestamp = clock();
	start_time = clock();
}

//prints how long the program has executed and time since previous timestamp
void printTimestamp()
{
	double total_time_elapsed  = (std::clock() - start_time) / (double) CLOCKS_PER_SEC;
	double time_since_previous = (std::clock() - timestamp) / (double) CLOCKS_PER_SEC;
		
	printf("TIME ELAPSED: %.02f\t(%.02f since previous)\n", total_time_elapsed, time_since_previous);
	
	timestamp = clock();
}

//prints a standardized header for the output csv
void createHeader(string benchmark_name)
{
	//print to output .csv file
	string output_file_name = "./results/" + benchmark_name + ".csv";
	FILE* out_file = fopen(output_file_name.c_str(), "w"); //open csv for writing
	
	fprintf(out_file, "\t\tWirelength\t\t\tArea\t\n");
	fprintf(out_file, "Trial\tnum_shapes\tInitial\tFinal\t%% change\tBest\tFinal\t %% from best\tExecution Time(s)\n");
				
	fclose(out_file);
}


//print the results of the full annealing run to the console and to a csv file
void printAnnealingResults(int trial_number, string benchmark_name, float initial_wirelength)
{
	float best_possible_area = bestPossibleArea();
	float final_area = best_layout->area;
	float area_percent = (float)(final_area - best_possible_area) / (float)best_possible_area;
	
	float final_wirelength = estimateTotalWirelength(best_layout);
	float wirelength_percent_change = (float)(initial_wirelength - final_wirelength) / (float)initial_wirelength;
	float time_elapsed = (std::clock() - start_time) / (float) CLOCKS_PER_SEC;
			
	//print to console
	printf("RUN #%d COMPLETE**********\nRESULTS FOR %s:\n", trial_number, benchmark_name.c_str());	
//printf("effort: %.2f\t", effort);
	printf("Wirelength:\t Initial \t Final \t %% change\n");
	printf("\t\t %8.0f \t %8.0f \t %.2f\n\n", initial_wirelength, final_wirelength, wirelength_percent_change);
	printf("Area:\t\t Best \t\t Final \t %% from best\n");
	printf("\t\t %8.0f \t %8.0f \t %.2f\n\n", best_possible_area, final_area, area_percent);
	printf("TIME ELAPSED: %.1f sec\n", time_elapsed);
	printf("\n");
	
	
	//print to output .csv file
	string output_file_name = "./results/" + benchmark_name + ".csv";
	FILE* out_file = fopen(output_file_name.c_str(), "a"); //open csv for appending
	
	
	fprintf(out_file, "%d\t%d\t%8.0f\t%8.0f\t%.2f\t%8.0f\t%8.0f\t%.2f\t%.1f\n", trial_number, num_shapes, initial_wirelength, final_wirelength, wirelength_percent_change, best_possible_area, final_area, area_percent, time_elapsed);
			
	fclose(out_file);

	resetTimestamp();
}

//print statistics for a specific BLOCK
void printBlockShapes(BLOCK* b)
{	
	printf("Block %s area = %.3f\tnum shapes carried: %d\n", b->name.c_str(), b->area, b->shapes.size());
	
	for(int j = 0; j < b->shapes.size(); j++)
	{
		printf("first: %.3f second: %.3f\n", b->shapes[j].first, b->shapes[j].second);	
	}
	
	printf("width: %.3f height: %.3f\n", b->width, b->height);
	printf("x: %.3f y: %.3f\n", b->x, b->y);
		
	printf("\n");
}

//returns the sum of all block areas
float bestPossibleArea()
{
	float sum = 0;
	for(int i = 0; i < blocks.size(); i++)
		sum += blocks[i]->area;
		
	return sum;
}

//export the block (usually the whole layout) to a text file which gives shapes for each block
void exportBlockShapes(BLOCK* b, string name)
{
	solidifyBlocks(b);
	
	FILE *output_file;
	string path = "./layouts/layout_" + name;

	output_file = fopen(path.c_str(), "w");	
	if (output_file == NULL) { fprintf(stderr, "Can't open output_layout file!");	exit(1); }
	
	printf("\nExporting layout file. Area = %.1f \t Wirelength = %.1f \t Combined Cost = %.1f \n\n", b->area, estimateTotalWirelength(b), combinedCost(b));
	
	fprintf(output_file, "name\tx\ty\tw\th\n");
		
	printBlockToFile(b, output_file);
	
	fclose(output_file);
}

void printBlockToFile(BLOCK* b, FILE* output_file)
{
	if(b->sub_blocks.size() == 0) //if no sub-blocks, it is basic block. Print to file
	{
		fprintf(output_file,"%s\t%.3f\t%.3f\t%.3f\t%.3f\n", b->name.c_str(), b->x, b->y, b->width, b->height);
		return;
	}
	else //else it has lower-level blocks. Print those.
	{
		//fprintf(output_file,"%s\t%.3f\t%.3f\t%.3f\t%.3f\n", b->name, b->x, b->y, b->width, b->height);
		for(int i = 0; i < b->sub_blocks.size(); i++)
			printBlockToFile(b->sub_blocks[i], output_file);
	}
}

float initial_wirelength;

void anneal_slicing(string benchmark_name, int trial_num)
{
	printf("\n******BEGIN ANNEALING RUN #%d******\n", trial_num);

	string benchmark_path = "./fpga_benchmarks/" + benchmark_name + ".fpga";
	
	//read the design files and store in vectors
	readBlocksFile(benchmark_path);				
	initializeOps(); //operands: H = Horizontal cut, V = Veritical cut
	
	initializeAnnealing();
		
	bestArea_time = 0;
		
	int move_count = 0;
	
	float temp = initializeTemp();
	float next_temp;
	
		
	
	initial_wirelength = estimateTotalWirelength(best_layout);
	
	float print_temp = temp;

	num_moves_accepted = 0;
	num_moves_accepted_Pr = 0;
	num_moves_rejected = 0;
	accepted_deltas = 0;
	accepted_Pr_deltas = 0;
	rejected_deltas = 0;
	
	printf("\nStarting layout: %d blocks\n", blocks.size());
	printf("Area: %.3f \t W: %.3f \t H: %.3f\t Wirelength: %.3f\n", best_layout->area, best_layout->width, best_layout->height, initial_wirelength);
	for(int i = 0; i < blocks.size(); i++)	
		printf("%s '%s' ", blocks[i]->name.c_str(), ops[i].c_str());	
	printf("\n");	
		
	while(temp > 0.1) //stop at temp = 0.1
	{						
		for(int i = 0; i < STEPS_PER_TEMP; i++)
			performAnnealingStep(temp);
		
		if(temp <= print_temp) //print the current state for every 90% reduction in temp	
		{	
			printf("########################\nTEMP = %.1f\talpha = %.2f\tbeta = %.2f\n", temp, alpha, beta);
			printf("num_active_blocks = %d\tbestArea_time = %.1f\n", num_active_blocks, bestArea_time/ (double) CLOCKS_PER_SEC);
			float average = accumulate( all_shapes_carried_forward.begin(), all_shapes_carried_forward.end(), 0.0)/all_shapes_carried_forward.size(); 
			printf("Avg # shapes carried forward: %.1f\n", average);
			all_shapes_carried_forward.clear(); //reset
			
			printf("Area: %.1f\tW: %5.1f\tH: %5.1f \t Wirelength: %.1f\n", best_layout->area, best_layout->width, best_layout->height, estimateTotalWirelength(best_layout));
			
			printf("# accepted: %d \t # accepted w/ Pr: %d \t # rejected: %d\n", num_moves_accepted, num_moves_accepted_Pr, num_moves_rejected);
			
			printf("avg accepted deltas: %.3f \t avg accepted_Pr deltas: %.3f \t avg rejected deltas: %.3f\n", (float)(accepted_deltas / (num_moves_accepted)), (float)(accepted_Pr_deltas / (num_moves_accepted_Pr)), (float)(rejected_deltas / (num_moves_rejected)) );
	
				
			num_moves_accepted = 0;
			num_moves_accepted_Pr = 0;
			num_moves_rejected = 0;
			accepted_deltas = 0;
			accepted_Pr_deltas = 0;
			rejected_deltas = 0;
			
			printTimestamp();
			print_temp *= .1;						
		}
		
		
		next_temp = temp * REDUCTION_FACTOR;
		
		//check for change in parameters based on temperature
		if(temp > 1 && next_temp <= 1)
		{
			//tip the weights towards area below 100
			//alpha = 0.9;
			//beta  = 0.1;
		}
		
		temp = next_temp;
	}
	
	//finished annealing, print results!

	printf("Final: %d blocks\n", blocks.size());
	printf("Area: %.3f \t W: %.3f \t H: %.3f\n", best_layout->area, best_layout->width, best_layout->height);
	for(int i = 0; i < blocks.size(); i++)	
		printf("%s '%s' ", blocks[i]->name.c_str(), ops[i].c_str());
		
	string name = benchmark_name + "_" + to_string((long long int) num_shapes) + "shapes_trial" + to_string((long long int) trial_num);
	
	exportBlockShapes(best_layout, name);
	
	//deallocate memory of the blocks for this run.
	//deallocateBlock(best_layout);
	
	return;
}


int main(int argc, char *argv[])
{	


//generate base blocks from the input file

	srand(time(0)); //use current time for random seed
	
	
	string benchmark_name = "hpnew"; //default
	if(argc > 1) benchmark_name = argv[1]; //arg1 is the name of the benchmark to use	
	
	if(argc > 2) num_shapes = stoi(argv[2]); //arg2 is the number of shapes to examine for each block
	
	//anneal_slicing(benchmark_path); printAnnealingResults(1, benchmark_name, initial_wirelength); return 1;
	
		
//run all benchmarks through many different number of shapes, 5 trials each
		
	vector <string> benchmark_names = {
	//"xeroxnew", 
	//"aptenew", 
	//"ami33new", 
	//"ami49new", 
	"n100anew", 
	"n200anew", 
	"n300anew",
	//"hpnew"
	};
	
	vector <int> num_shapes_to_try = {
	//5,
	//9,
	//15,
	19,
	//25
	
	};
	
	int num_runs = 3;
	
	for(int benchmark_index = 0; benchmark_index < benchmark_names.size(); benchmark_index++)
	{
		benchmark_name = benchmark_names[benchmark_index];		
		
		createHeader(benchmark_name); //insert the header into the results csv
		
		for(int num_shapes_index = 0; num_shapes_index < num_shapes_to_try.size(); num_shapes_index++)
		{
			num_shapes = num_shapes_to_try[num_shapes_index];
							
			//loop through multiple annealing runs
			for(int trial_num = 1; trial_num <= num_runs; trial_num++)
			{							
				anneal_slicing(benchmark_name, trial_num);	
				printAnnealingResults(trial_num, benchmark_name, initial_wirelength);
			}
		}
	}
} //end main
