//anneal_std_cell.cpp
//Mark Sears and Navya Sreeram
//PA5: Layout of fixed standard cells using Simulated Annealing process

#include <stdio.h>
#include <stdlib.h>  
#include <stdint.h>
#include <string.h>
#include <string>
#include <list>
#include <map> 
#include <vector> 
#include <iterator>
#include <iostream>
#include <fstream>
#include <ctime>
#include <math.h>

using namespace std;

class NET;
class CELL;
class PAD;

//***************************//
//  DECLARE GLOBAL VARIABLES //
//***************************//

string file_path = "./IBM_Benchmarks/"; //path where design files can be found

vector <CELL*> cells; //pointers to all cells in the design

map <string, CELL*> cells_map; //map of cells. Key = cell_name, value = cell pointer
				//useful for accessing cells by name
map <string, PAD*> pads_map;
vector <NET*> nets; //pointers to all nets in the design

int CELL_HEIGHT = 16;
int CELL_WIDTH = 4;

double ROW_SPACING = 16; //inter-row die spacing fixed at 16
double PAD_SPACING = 4; //allowed spacing between 2 pads

double DIE_WIDTH;
double DIE_HEIGHT;
int NUM_ROWS; //number of rows in the layout
int ROW_LENGTH; //physical size of a row


std::clock_t start_time;
std::clock_t timestamp;


//***************************//
//     FUNCTION HEADERS	     //
//***************************//

class NET {
	public:
	char name[100];
	int num_pins; //total number of pins on this net
	bool cutstate; //true if the net is cut across the partition
	bool critical; //a net if critical if there is a cell move that affects cutstate
	bool has_locked[2];//a net is "dead" if no legal cell swaps remain that could change the cutstate
			//this occurs when there is a locked or fixed cell in both partitions on this net
			
	int partition_count[2]; //count of how many cells in each partition on this net
	
	vector <CELL*> cell_list; //pointers to all cells on this net
	vector <PAD*> pad_list; //pointers to all pads on this net			
};

class CELL {
	public:
	char name[100];
	char type[20];
	
	int width;
	int height;
	int area; 	//physical area of the cell
	
	int x;	//(x,y) coordinates of the cell's lower left corner
	int y;
	
	/*
	bool fixed; 	//true for cells that can never be moved
	bool locked; 	//true for cells locked for the rest of the pass
	int gain; 	//gain is the change in cutset if this cell were to swap partitions
	int partition; 	//which partition this cell is in: 0 or 1
	int best_partition; //used to remember this cell's partition in the best solution found so far
	multimap <int, CELL*> :: iterator gain_itr; //an iterator that points to this cell's location in the gain buckets	
	*/
	
	vector <NET*> net_list; //pointers to all nets this cell is on
			
};

class PAD { 	//pads have no size but are used to denote a wire connection
	public:
	char name[100];
	int x;	//(x,y) coordinates of the cell's lower left corner
	int y; 
	
	vector <NET*> net_list; //pointers to all nets this pad is on
};

//reads the nets file "design.nets". Fills in "nets" array with all the nets in the design.
void readNetsFile(string nets_file_name)
{	
	nets.clear();
	cells.clear();
	pads_map.clear();
	cells_map.clear();
	
	ifstream nets_file;			
	
	nets_file.open(nets_file_name.c_str());			
	if (nets_file == NULL) { fprintf(stderr, "Can't open nets file: %s\n", nets_file_name.c_str());	exit(1); }
	
	printf("\nReading nets file: %s\n", nets_file_name.c_str());	
	
	char next[20]; //reads the next word of characters
	char cell_name[20];
		
	for(int i = 0; i < 5; i++)	//ignore the first 5 lines
		nets_file >> next;					
	
	NET* current_net = new NET;
	
	while(1)
	{
		nets_file >> cell_name;
				
		if(nets_file.eof())
		{			
			nets.push_back(current_net);
			break;
		}
			
		if(cell_name[0] == 'p' || cell_name[0] == 'a') //read in a cell (this should always trip)
		{												
			if(cell_name[0] == 'a')
			{
				if(cells_map[cell_name] == NULL) //if this cell doesn't exist, create it
				{
					CELL* current_cell = new CELL;
					current_cell->width  = 4; //fixed cell dimensions
					current_cell->height = CELL_HEIGHT;					
					strcpy(current_cell->name, cell_name);					
					cells_map[cell_name] = current_cell;
					cells.push_back(current_cell);
				}
			}
					
			if(cell_name[0] == 'p')
			{				
				if(pads_map[cell_name] == NULL) //if this pad doesn't exist, create it
				{	
					PAD* current_pad = new PAD;						
					strcpy(current_pad->name, cell_name);					
					pads_map[cell_name] = current_pad;
				}
			}
									
			nets_file >> next;
			if (next[0] == 's') //at the start of new nets
			{				
				if(current_net->cell_list.size() != 0) //if this net has cells
					nets.push_back(current_net); //add the previous net to the set
					
				current_net = new NET; 				
				
				nets_file >> next; //read in the extra character
			}
			
			if(cell_name[0] == 'a')
			{
				current_net->cell_list.push_back(cells_map[cell_name]); //add this cell to the net
				cells_map[cell_name]->net_list.push_back(current_net);
			}
			else if(cell_name[0] == 'p')
			{
				current_net->pad_list.push_back(pads_map[cell_name]); //add this pad to the net
				pads_map[cell_name]->net_list.push_back(current_net);
			}
			
		} else printf("ERROR reading nets!\n");		 
	}	
	
	nets_file.close();
	printf("Done!\n\n");
	return;			
}

//initializePlacement() initilizes the placement of all cells by setting coordinates for lower-left corner
//such that the layout has an aspect ratio = 1
void initializePlacement()
{
	int total_cell_width = 0; //add up all the cell widths
	for (std::map<string, CELL*>::iterator it=cells_map.begin(); it!=cells_map.end(); ++it)
		total_cell_width += it->second->width;		

	//initialize chip dimensions, aim for Aspect Ratio = 1
	NUM_ROWS = ceil(sqrt(((double)total_cell_width) / ((double)(CELL_HEIGHT + ROW_SPACING))));
	ROW_LENGTH = (int)((double)total_cell_width / (double)NUM_ROWS);	
	ROW_LENGTH -= (ROW_LENGTH%4); //make sure ROW_LENGTH is a multiple of 4
	
	DIE_WIDTH = ROW_LENGTH;
	DIE_HEIGHT = NUM_ROWS * (CELL_HEIGHT + ROW_SPACING);	
	
	//try to set PAD_SPACING such that all the pads fit on 4 sides
	PAD_SPACING = floor((double)(4 * DIE_WIDTH) / (double)(pads_map.size() + 4));
	printf("PAD_SPACING: %.2f\n", PAD_SPACING);
	
	bool sides_full = false; //set true after left and right sides are full up on pads	
	double next_x = 0; //coordinates for placing cells
	double next_y = PAD_SPACING; //don't let pads sit in the corners
	
	for (std::map<string, PAD*>::iterator it=pads_map.begin(); it!=pads_map.end(); ++it)
	{
		PAD* current_pad = it->second;
		
		current_pad->x = next_x;
		current_pad->y = next_y;
		
		if(!sides_full) //place pads on alternating right and left sides
		{
			if(next_x == ROW_LENGTH)
			{
				next_x = 0;
				next_y += PAD_SPACING;
			} else next_x = ROW_LENGTH;
			
			if(next_y >= DIE_HEIGHT)
			{
				sides_full = true;
				next_y = 0;
				next_x = PAD_SPACING;
			}
		}
		else	//sides full, place pads on alternating top and bottom sides
		{
			if(next_y == DIE_HEIGHT)
			{
				next_y = 0;
				next_x += PAD_SPACING;
			} else next_y = DIE_HEIGHT;
			
			if(next_x >= ROW_LENGTH)
			{
				printf("ERROR during initial placement! Too many pads!\n");
			}
		}
	}	
	
	next_x = 0; //reset
	next_y = 0;
	
	//for each cell, give it a placement within a row
	for (std::map<string, CELL*>::iterator it=cells_map.begin(); it!=cells_map.end(); ++it)
	{
		CELL* current_cell = it->second;
		
		current_cell->x = next_x;
		current_cell->y = next_y;
		
		next_x += CELL_WIDTH;
		
		if(next_x > ROW_LENGTH)
		{
			next_x = 0;
			next_y += CELL_HEIGHT + ROW_SPACING;
			
			if(next_y > NUM_ROWS * (CELL_HEIGHT + ROW_SPACING)) 
				printf("ERROR during initial placement! Too many cells!\n");					
		}
	}
	
}



//returns the half perimeter or "bounding box" wire length esimate for a net
//includes pads 
double halfPerimeter(NET* net)
{
	double min_x = net->cell_list[0]->x, max_x = net->cell_list[0]->x;
	double min_y = net->cell_list[0]->y, max_y = net->cell_list[0]->y;
	
	//find the min and max of the bounding box
	for(int i = 1; i < net->cell_list.size(); i++)
	{
		if(net->cell_list[i]->x < min_x)
			min_x = net->cell_list[i]->x;
		if(net->cell_list[i]->x > max_x)
			max_x = net->cell_list[i]->x;
		
		if(net->cell_list[i]->y < min_y)
			min_y = net->cell_list[i]->y;
		if(net->cell_list[i]->y > max_y)
			max_y = net->cell_list[i]->y;
	}
	
	//include the pads!
	for(int i = 0; i < net->pad_list.size(); i++)
	{
		if(net->pad_list[i]->x < min_x)
			min_x = net->pad_list[i]->x;
		if(net->pad_list[i]->x > max_x)
			max_x = net->pad_list[i]->x;
		
		if(net->pad_list[i]->y < min_y)
			min_y = net->pad_list[i]->y;
		if(net->pad_list[i]->y > max_y)
			max_y = net->pad_list[i]->y;
	}
	
	return (max_x - min_x) + (max_y - min_y);
}


void swapCellCoords(CELL* c1, CELL* c2)
{
	int temp_x = c1->x;
	int temp_y = c1->y;
	
	c1->x = c2->x;
	c1->y = c2->y;
	
	c2->x = temp_x;
	c2->y = temp_y;
}

//swaps the physical position of 2 given cells
//returns the change in cost of performing the move
double swapCellPlacement(CELL* c1, CELL* c2)
{
	double initial_cost = 0; //cost before the move
	for(int i = 0; i < c1->net_list.size(); i++)
		initial_cost += halfPerimeter(c1->net_list[i]);
	for(int i = 0; i < c2->net_list.size(); i++)
		initial_cost += halfPerimeter(c2->net_list[i]);
	
	swapCellCoords(c1, c2); //perform the coordinate swap
	
	double final_cost = 0; //cost after the move
	for(int i = 0; i < c1->net_list.size(); i++)
		final_cost += halfPerimeter(c1->net_list[i]);
	for(int i = 0; i < c2->net_list.size(); i++)
		final_cost += halfPerimeter(c2->net_list[i]);

		
	return final_cost - initial_cost;
}


//returns the total estimated wire length for the whole layout
double findTotalCost()
{
	double total_wirelength = 0;
	//for each net, estimate wirelength
	for(int i = 0; i < nets.size(); i++)
	{
		total_wirelength += halfPerimeter(nets[i]);
	}
	return total_wirelength;
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

//returns an intial value for temperature
//perform a "dry run" of many annealing steps to determine the average change
double initializeTemp()
{
	vector <double> changes;
	int num_swaps = 10000;
	
	//perform random cell swaps
	for(int i = 0; i < num_swaps; i++)	
		changes.push_back(swapCellPlacement(cells[rand()%cells.size()], cells[rand()%cells.size()]));	
	
	double std_dev = StandardDeviation(changes);
		
	return 30*std_dev; //for a temperature of 30*std_dev,
			//there will be a 90% chance to accept a very bad change of 3*std_dev
}

double performAnnealingStep(double temp)
{
	CELL* c1 = cells[rand()%cells.size()]; //selected 2 random cells to swap
	CELL* c2 = cells[rand()%cells.size()];
	double delta = swapCellPlacement(c1, c2);
		
	if(delta <= 0) //if the change is good, accept it
		return delta;
	
	//otherwise, accept with probability		
	double Pr_accepting = exp((double)(-1*delta)/(double)temp);
	int p = 100000;
	if((double)(rand()%p) / (double)p > Pr_accepting) //based on temp and change, chance to accept
	{
		swapCellCoords(c1, c2); //reject the change, swap the cells back
		delta = 0;
	}
	//otherwise, the change is accepted.
	
	return delta;
}

//prints how long the program has executed and time since previous timestamp
void printTimestamp()
{

	printf("TIME ELAPSED: %.02f\t(%.02f since previous)\n", (std::clock() - start_time) / (double) CLOCKS_PER_SEC, (std::clock() - timestamp) / (double) CLOCKS_PER_SEC);
	
	timestamp = clock();
}

void resetTimestamp()
{
	timestamp = clock();
	start_time = clock();
}

void printDimensions()
{
	printf("Num_cells: %d\tNum_pads: %d\tNum_nets: %d\n", cells.size(), pads_map.size(), nets.size());
	printf("DIE_WIDTH: %.2f\tDIE_HEIGHT: %.2f\tASPECT RATIO (W/H): %.2f\n", DIE_WIDTH, DIE_HEIGHT, DIE_WIDTH/DIE_HEIGHT);
	printf("NUM_ROWS: %d\tROW_LENGTH: %d\n", NUM_ROWS, ROW_LENGTH);
}

void printCellInfo()
{
	printf("CELLS:\n");for (std::map<string, CELL*>::iterator it=cells_map.begin(); it!=cells_map.end(); ++it)    	
			printf("%s \t(%d,%d)\t is on %d nets.\n", it->first.c_str(), it->second->x, it->second->y, it->second->net_list.size());
					
	
	printf("PADS:\n"); for (std::map<string, PAD*>::iterator it=pads_map.begin(); it!=pads_map.end(); ++it)
		printf("%s \t(%d,%d)\t is on %d nets.\n", it->first.c_str(), it->second->x, it->second->y, it->second->net_list.size());
	
		
	printf("\nNETS:\n");
	for(int i = 0; i < nets.size(); i++)
    	{
		printf("Net %d has cells:", i);
		for(int j = 0; j < nets[i]->cell_list.size(); j++)
			printf(" %s", nets[i]->cell_list[j]->name);
		for(int j = 0; j < nets[i]->pad_list.size(); j++)
			printf(" %s", nets[i]->pad_list[j]->name);
			
		printf("\tHalf Perimeter: %.2f", halfPerimeter(nets[i]));
		printf("\n");	
	}
}

double performAnnealing(int run_number, string benchmark_name, double effort)
{
	printf("\n**********STANDARD CELL PLACEMENT ANNEALING**********\n");
	printf("**********RUN #%d**********\n", run_number);
	
	string nets_file_name = file_path + benchmark_name + "/" + benchmark_name + ".net";	
	readNetsFile(nets_file_name);
	initializePlacement();
	printDimensions();
	
	double temp = initializeTemp();		
	double initial_wirelength = findTotalCost();
	double current_wirelength = initial_wirelength;	
	
	int swaps_per_temp = nets.size() * 4 * log10((effort+1)); //set the number of swaps at each temp
	double REDUCTION_FACTOR = 1 - .1*exp(-0.2*effort);
	
	printf("effort: %.2f\tswaps_per_temp: %d\tREDUCTION_FACTOR: %.3f\n", effort, swaps_per_temp, REDUCTION_FACTOR);
	if(effort > 10) printf("###### WARNING! High effort may cause excessively long execution time!\n");
	
	
	bool frozen = false; //annealing end condition
	double print_temp = temp;
		
	printf("\n**********BEGIN ANNEALING**********\n");	
	do {
		if(temp <= print_temp) //print the current state for every 90% reduction in temp	
		{
			printf("Wirelength: %.0f\tTemp: %.2f\t", current_wirelength, temp);
			printTimestamp();
			print_temp *= .1;
		}
		
		int num_swaps_accepted = 0;
		
		//perform swaps, then decrease temp
		for(int i = 0; i < swaps_per_temp; i++)
		{
			double delta = performAnnealingStep(temp);
			current_wirelength += delta;
			
			if(delta != 0) //count the number of accepted steps
				num_swaps_accepted++;
		}
				
		
		temp *= REDUCTION_FACTOR;
		
		//stop at low temp or if less than 0.1% of swaps are accepted		
		if(temp < .1)
		{			
			printf("\n**********STOPPING due to low temp: %.2f**********\n", temp);
			frozen = true;
		}		
		if(num_swaps_accepted * 1000 < swaps_per_temp) 
		{
			printf("\n**********STOPPING due to < 0.1%%  swaps accepted: %d\n", num_swaps_accepted);
			frozen = true;
		}
	} while(!frozen);	
	
	double final_wirelength = findTotalCost();
	double percent_change = 100*(initial_wirelength - final_wirelength) / initial_wirelength;
	double time_elapsed = (std::clock() - start_time) / (double) CLOCKS_PER_SEC;	
			
	//print to console
	printf("RUN #%d***RESULTS FOR %s:\n", run_number, (file_path + benchmark_name + "/").c_str());	
	printf("effort: %.2f\t", effort);
	printf("Final Wirelength: %.0f\t", final_wirelength);
	printf("%% change: %.1f\t", 100*(initial_wirelength - final_wirelength) / initial_wirelength );
	printf("TIME ELAPSED: %.1f sec\n", time_elapsed);
	printf("\n");
	
		
	//print to output .csv file
	string output_file_name = file_path + benchmark_name + "/anneal_output.csv";
	FILE* out_file = fopen(output_file_name.c_str(), "a");
	
	
	fprintf(out_file, "%d,\t%.2f,\t%.0f,\t%.1f,\t\t%.1f\n", run_number, effort, final_wirelength, percent_change, time_elapsed);
			
	fclose(out_file);
	
	resetTimestamp();
	return final_wirelength;
}

int main(int argc, char *argv[])
{	
	start_time = std::clock();
	timestamp = std::clock();
	srand(time(0)); //use current time for random seed
	
	//read the example number from the argument to the filepath
	string benchmark_name = "ibm";	
	if(argc > 1) benchmark_name += argv[1];
	else benchmark_name += "01"; //defaults to example 01	
	
	//read the effort parameter
	double effort;	
	if(argc > 2) effort = stof(argv[2]);
	else effort = 6.0; //default
				
	int num_runs = 1;	
	for(int run = 1; run <= num_runs; run++)			
		performAnnealing(run, benchmark_name, effort);	
	/*
	//perform multiple sample runs on all benchmarks
	//record results in text files
	int num_runs = 1;
	string benchmark_name = "ibm";	
	string benchmark_nums[10] = {"01", "03", "04", "06", "08", "10", "12", "14", "16", "18"};
	//string benchmark_nums[4] = {"01", "03", "08", "10"};
	
	//"effort" quantifies how hard the algorithm will work to find a good solution
	//effort can be any positive number, but is designed to give best results when  1 <= effort <= 10
	//effort = 1 gives poor total wirelength, but fast execution time
	//effort = 10 gives good total wirelength but takes longer
	//efforts above 10 do not see significant improvement to wirelength, but still increases execution time	
	for(int a = 0; a < 1; a++)
	{
		string output_file_name = file_path + benchmark_name + benchmark_nums[a] + "/anneal_output.csv";
			FILE* out_file= fopen(output_file_name.c_str(), "w"); //print column headers		 
			fprintf(out_file, "Benchmark: %s\n", (benchmark_name + benchmark_nums[a]).c_str());
			fprintf(out_file, "Run,\tEffort,\tFinal Cost,\t%% change,\tTime (s)\n");
			fclose(out_file);
			
		//for(int i = 0; i < 6; i++)
		for(double effort = 7; effort <= 10; effort += .5)		
			for(int run = 1; run <= num_runs; run++)						
				performAnnealing(run, benchmark_name + benchmark_nums[a], effort);
	}
	*/				
}
