#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "mc.h"

using namespace std;

void mc :: read_from_in(ifstream& in)
{
	string label_max_iter = "max_iter";
	string label_T = "temperature";
	string label_if_test = "if_test";
	string label_act_p = "begin_action_probability";
	string label_num_ele = "num_ele";
	string tmp;
	stringstream ss;

	// get max number of iteration
	getline(in, tmp);
	while(tmp.find(label_max_iter) == string::npos)
		getline(in,tmp);
	ss << (tmp);
	getline(ss,tmp,'='); ss >> max_iter;
	ss.str(""); ss.clear(); in.clear(); in.seekg(ios::beg);

	// get running temperature
	getline(in, tmp);
	while(tmp.find(label_T) == string::npos)
		getline(in,tmp);
	ss << (tmp);
	getline(ss,tmp,'='); ss >> T;
	ss.str(""); ss.clear(); in.clear(); in.seekg(ios::beg);

	// get if running test mode (No QE call)
	getline(in, tmp);
	while(tmp.find(label_if_test) == string::npos)
		getline(in,tmp);
	ss << (tmp);
	getline(ss,tmp,'='); ss >> if_test;
	ss.str(""); ss.clear(); in.clear(); in.seekg(ios::beg);

	// get action probability
	getline(in, tmp);
	while(tmp.find(label_act_p) == string::npos)
		getline(in,tmp);
	getline(in, tmp);
	ss << (tmp);
	for(size_t t1=0; t1<num_act; t1++)
		act_p[t1]=0;
	for(size_t t1=0; t1<num_act; t1++)
		ss>>act_p[t1];
	ss.str(""); ss.clear(); in.clear(); in.seekg(ios::beg);

	// get number of elements
	getline(in,tmp);
	while(tmp.find(label_num_ele) == string::npos)
		getline(in,tmp);
	ss << (tmp);
	getline(ss,tmp,'='); ss >> num_ele;
	ss.str(""); ss.clear(); in.clear(); in.seekg(ios::beg);

	// initialize other parameters
	num_atm_each_change.resize(num_ele);
	act_type = -1;
	opt_e = 0;
}

void mc :: create_new_structure(cell c_old, cell& c_new)
{
	c_new = c_old;
	double tmp_p[3];
	// for examine add
	int num_trial = 10000;
	vec pos_add;
	int ele_type_add;
	int atm_id_tmp, atm_id_tmp2;
	atom atm_tmp;
	element ele_tmp;
	double r_tmp;
	double add_p_sum;
	// for exam remove and swap
	int num_removable_ele;
	int num_movable_ele;
	//=======================================
	// check if each action is accessable, and adjust p if not
	cout<<"Begin adjust weight of actions:"<<endl;
	//---------------------------------------
	// check for add
	if (act_p[0] > 0)
	{
		// choose type of element to add
		add_p_sum = 0;
		for(size_t t1=0; t1<c_new.num_ele; t1++)
			add_p_sum += c_new.ele_list[t1].p_add;
		add_p_sum = (double)rand()/RAND_MAX * add_p_sum;
		for(size_t t1=0; t1<c_new.num_ele; t1++)
		{
			if( add_p_sum <= c_new.ele_list[t1].p_add)
			{
				ele_type_add = t1;
				ele_tmp = c_new.ele_list[t1];
				break;
			}
			else
			{
				add_p_sum -= c_new.ele_list[t1].p_add;
			}
		}
		for(size_t t1=0; t1<num_trial; t1++)
		{
			// choose where attemp to add
			atm_id_tmp = rand()%c_new.num_atm;
			atm_tmp = c_new.atm_list[atm_id_tmp];
			pos_add = atm_tmp.pos + pos_add.rand_norm()*(ele_tmp.r_min + (double)rand()/RAND_MAX*(ele_tmp.r_max - ele_tmp.r_min));
			c_new.min_distance(pos_add, r_tmp, atm_id_tmp);
			// if do vc-relax, then only need to satisfy coordination rule, otherwise new position has to be between r_min and r_max
			if (r_tmp>ele_tmp.r_min && r_tmp<ele_tmp.r_max && (c_new.if_vc_relax || (pos_add.x[2] > c_new.h_min && pos_add.x[2] < c_new.h_max)))
				break;
		}
		if (r_tmp>ele_tmp.r_min && r_tmp<ele_tmp.r_max && (c_new.if_vc_relax || (pos_add.x[2] > c_new.h_min && pos_add.x[2] < c_new.h_max)))
			tmp_p[0] = act_p[0];
		else
		{
			tmp_p[0] = 0;
			cout<<"    Can not find site to add atom, weight of choosing to add is set to 0"<<endl;
		}
	}
	else
		tmp_p[0] = 0;
	//---------------------------------------
	// check for remove
	if (act_p[1] > 0)
	{
		num_removable_ele = 0;
		for (size_t t1=0; t1<c_new.num_ele; t1++)
		{
			if (c_new.num_ele_each_remove[t1] > 0)
				num_removable_ele++;
		}
		if (num_removable_ele >= 1)
			tmp_p[1] = act_p[1];
		else
		{
			tmp_p[1] = 0;
			cout<<"    Can not find removable elements, weight of choosing to remove set to 0"<<endl;
		}
	}
	else
		tmp_p[1] = 0;
	// check for swap
	if (act_p[2] > 0)
	{
		num_movable_ele = 0;
		for (size_t t1=0; t1<c_new.num_ele; t1++)
		{
			if (c_new.num_ele_each_move[t1] > 0)
				num_movable_ele++;
		}
		if (num_movable_ele >= 2)
			tmp_p[2] = act_p[2];
		else
		{
			tmp_p[2] = 0;
			cout<<"    Can not find more than one movable elements, weight of choosing to swap is set to 0"<<endl;
		}
	}
	else
		tmp_p[2] = 0;
	cout<<"    New weight of each action is:"<<endl;
	cout<<"    Add: "<<tmp_p[0]<<"    Remove: "<<tmp_p[1]<<"    Swap: "<<tmp_p[2]<<endl;
	cout<<"End adjust weight of actions"<<endl<<endl;
	//=======================================
	// create new structure
	// decide which action to choose
	add_p_sum = 0;
	for(size_t t1=0; t1<3; t1++)
		add_p_sum += tmp_p[t1];
	// exit if no actions are allowed
	if (fabs(add_p_sum) <= 1e-10)
	{
		cout<<"Error: No actions are legal"<<endl;
		exit(EXIT_FAILURE);
	}
	add_p_sum = (double)rand()/RAND_MAX * add_p_sum;
	act_type = -1;
	for(size_t t1=0; t1<3; t1++)
	{
		if(add_p_sum < tmp_p[t1])
		{
			act_type=t1;
			break;
		}
		else
		{
			add_p_sum -= tmp_p[t1];
		}
	}
	// start applying change
	//=======================================
	for (size_t t1=0; t1<num_ele; t1++)
		num_atm_each_change[t1] = 0;
	switch(act_type)
	{
		//---------------------------------------
		// add
		case 0:
		{
			num_atm_each_change[ele_type_add] = 1;
			// adjust pos_add so it is in the cell
			vec bb[3];
			double coef_bb[3];
			bb[0] = (c_new.latt[1]^c_new.latt[2])/((c_new.latt[1]^c_new.latt[2])*c_new.latt[0]);
			bb[1] = (c_new.latt[2]^c_new.latt[0])/((c_new.latt[1]^c_new.latt[2])*c_new.latt[0]);
			bb[2] = (c_new.latt[0]^c_new.latt[1])/((c_new.latt[1]^c_new.latt[2])*c_new.latt[0]);
			for(size_t t1=0; t1<3; t1++)
			{
				coef_bb[t1] = pos_add*bb[t1];
				coef_bb[t1] -= floor(coef_bb[t1]);
			}
			pos_add = c_new.latt[0]*coef_bb[0]+c_new.latt[1]*coef_bb[1]+c_new.latt[2]*coef_bb[2];
			c_new.ad_atom(pos_add,ele_type_add);
			cout<<"Atom added, +"<<c_new.num_atm<<" "<<c_new.ele_list[ele_type_add].sym<<", position is "<<pos_add<<endl;
			break;
		}
		//---------------------------------------
		// remove
		case 1:
		{
			// find which removable atom to remove
			atm_id_tmp = rand()%c_new.num_atm_remove;
			for(size_t t1=0; t1<c_new.num_atm; t1++)
			{
				if(atm_id_tmp == 0 && c_new.atm_list[t1].if_move == 2)
				{
					atm_id_tmp = t1;
					break;
				}
				else if (c_new.atm_list[t1].if_move == 2)
				{
					atm_id_tmp--;
				}
			}
			num_atm_each_change[c_old.atm_list[atm_id_tmp].type] = -1;
			c_new.rm_atom(atm_id_tmp);
			cout<<"Atom removed, -"<<atm_id_tmp+1<<" "<<c_old.atm_list[atm_id_tmp].ele->sym<<endl;
			break;
		}
		//---------------------------------------
		// swap
		case 2:
		{
			int iter_swap = rand()%c_old.num_atm+1;
			cout<<"Perform swap "<<iter_swap<<" times:"<<endl;
			for(size_t t_swap=0; t_swap<iter_swap; t_swap++)
			{
				// find the first atom to switch
				atm_id_tmp = rand()%c_new.num_atm_move;
				for(size_t t1=0; t1<c_new.num_atm; t1++)
				{
					if(atm_id_tmp == 0 && c_new.atm_list[t1].if_move >= 1)
					{
						atm_id_tmp = t1;
						break;
					}
					else if (c_new.atm_list[t1].if_move >= 1)
					{
						atm_id_tmp--;
					}
				}
				// find the second atom to switch
				atm_id_tmp2 = rand()%(c_new.num_atm_move - c_new.num_ele_each_move[c_new.atm_list[atm_id_tmp].type]);
				for(size_t t1=0; t1<c_new.num_atm; t1++)
				{
					if(atm_id_tmp2 == 0 && c_new.atm_list[t1].if_move >= 1 && c_new.atm_list[t1].type != c_new.atm_list[atm_id_tmp].type)
					{
						atm_id_tmp2 = t1;
						break;
					}
					else if (c_new.atm_list[t1].if_move >= 1 && c_new.atm_list[t1].type != c_new.atm_list[atm_id_tmp].type)
					{
						atm_id_tmp2--;
					}
				}
				c_new.sp_atom(atm_id_tmp,atm_id_tmp2);
				cout<<"Atoms swapped, "<<atm_id_tmp+1<<" "<<c_old.atm_list[atm_id_tmp].ele->sym<<" <--> "<<atm_id_tmp2+1<<" "<<c_old.atm_list[atm_id_tmp2].ele->sym<<endl;
			}
			break;
		}
		//---------------------------------------
		default:
		{
			cout<<"Error: Undifined action type: "<<act_type<<endl;
			exit(EXIT_FAILURE);
		}
	}
}

void mc :: save_opt_structure(cell& c_new)
{
	opt_e = c_new.energy;
	for(size_t t1=0; t1<c_new.num_atm; t1++)
		opt_e -= c_new.atm_list[t1].ele->mu;
	opt_c = c_new;
	cout<<"Initialized the minimum seeker to the starting structure"<<endl;
}

int mc :: check_if_accept(cell& c_old, cell& c_new)
{
	double exp_pre, exp_main;

	// calculate formation energy
	e1 = c_old.energy;
	e2 = c_new.energy;
	for(size_t t1=0; t1<c_old.num_atm; t1++)
		e1 -= c_old.atm_list[t1].ele->mu;
	for(size_t t1=0; t1<c_new.num_atm; t1++)
		e2 -= c_new.atm_list[t1].ele->mu;
	//==============================================
	cout<<endl<<"Evaluating whether to accept new structure"<<endl;
	if (fabs(c_new.energy) < 1e-9)
	{
		accept = 0;
		cout<<"    Warning: SCF of new structure does not converge"<<endl;
		cout<<"    ====Rejected===="<<endl;
	}
	else
	{
		switch (act_type)
		{
			//---------------------------------------
			// add,remove, and swap
			case 0:
			case 1:
			case 2:
			{
				// calculate prefactor and exp
				exp_pre = 1;
				c_new.get_volume();
				for(size_t t1=0; t1<num_ele; t1++)
				{
//					if(c_new.if_change_v)
//						exp_pre *= ( pow(c_new.vol,num_atm_each_change[t1])*factor(c_old.num_ele_each[t1])/pow(c_old.ele_list[t1].rho,-num_atm_each_change[t1])/factor(num_atm_each_change[t1]+c_old.num_ele_each[t1]) );
//					else
//					{
						c_old.ele_list[t1].update_tb(T);
						c_new.ele_list[t1].update_tb(T);
						exp_pre *= ( pow(c_new.vol,num_atm_each_change[t1])*factor(c_old.num_ele_each[t1])/pow(c_old.ele_list[t1].tb,3*num_atm_each_change[t1])/factor(num_atm_each_change[t1]+c_old.num_ele_each[t1]) );
//					}
				}
				// if allow V to change, multiply exp_pre by (V_new/V_old)^N_tot(old)
				if(c_new.if_change_v)
					exp_pre *= pow((c_new.vol/c_old.vol),c_old.num_atm);
				exp_main = exp((e1-e2)/T/kb);
				cout<<"    Pre. factor: "<<exp_pre<<"    exp. factor: "<<exp_main<<"    total factor: "<<exp_pre*exp_main<<endl;
				// start evaluating whether to accept or reject
				if (1 < exp_pre*exp_main)
				{
					accept = 1;
					cout<<"    ====Accepted===="<<endl;
					cout<<"    New structure has lower volume-factored-in formation energy"<<endl;
					if (e2 < opt_e)
					{
						opt_e = e2;
						opt_c = c_new;
						cout<<"    New structure has by far the lowest formation energy, best structure updated"<<endl;
					}
				}
				else if ((double)rand()/RAND_MAX < exp_pre*exp_main)
				{
					accept = 1;
					cout<<"    ----Accepted----"<<endl;
					cout<<"    New structure has highter volume-factored-in formation energy"<<endl;
					if (e2 < opt_e)
					{
						opt_e = e2;
						opt_c = c_new;
						cout<<"    New structure has by far the lowest formation energy, best structure updated"<<endl;
					}
				}
				else
				{
					accept = 0;
					cout<<"    ----Rejected----"<<endl;
				}
				break;
			}
			//---------------------------------------
			default:
			{
				cout<<"Error: Invalid action number: "<<act_type<<endl;
				exit(EXIT_FAILURE);
			}
		}
	}
	cout<<"Best formation energy is "<<setw(20)<<setprecision(9)<<opt_e<<" eV"<<endl;
	return accept;
}

int mc :: factor(int n)
{
	if (n < 0)
	{
		cout<<"Error: Number of atoms should not be negative"<<endl;
		exit(EXIT_FAILURE);
	}
	if (n==0)
		return 1;
	else
		return factor(n-1);
}

void mc :: print()
{
	cout<<"Max iteration: "<<max_iter<<endl;
	cout<<"Simulation temperature: "<<T<<endl;
	cout<<"If running test: "<<if_test<<endl;
	cout<<"Action probability: "<<endl;
	for(size_t t1=0; t1<num_act; t1++)
		cout<<'\t'<<act_p[t1];
	cout<<endl;
}
