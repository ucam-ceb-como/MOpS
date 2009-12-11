/*!
 * \file   swp_PAH_primary.h
 * \author Markus Sander
 *
 * \brief  Particle Model that stores the individual PAHs of a soot particle
 */
/*
  Author(s):      Markus Sander (ms785)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Markus Sander.

  File purpose:
    Implementation of the PAHPrimary class declared in the
    swp_PAH_primary.h header file.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Dr Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
*/

#include "swp_primary.h"
#include "swp_PAH_primary.h"
#include "swp_submodel_type.h"
#include "swp_aggmodel_type.h"
#include "swp_model_factory.h"
#include "swp_PAH_cache.h"
//#include <math.h>
#include "rng.h"
#include "swp_particle_image.h"
#include <stdexcept>
#include "string_functions.h"
using namespace Sweep;
using namespace Sweep::AggModels;
using namespace Sweep::SubModels;
using namespace std;
using namespace Strings;
/*
double PAHPrimary::pow(double a, double b) {   
    int tmp = (*(1 + (int *)&a));   
    int tmp2 = (int)(b * (tmp - 1072632447) + 1072632447);   
    double p = 0.0;   
    *(1 + (int * )&p) = tmp2;   
    return p;   
}  
*/

// CONSTRUCTORS AND DESTRUCTORS.
PAHPrimary::PAHPrimary() : Primary(),
    m_numcarbon(0),
    m_PAHmass(0),
    m_PAHCollDiameter(0),
    m_numPAH(0),
    m_numprimary(0),
    m_primarydiam(0.0),
    m_children_radius(0),
    m_children_vol(0),
    m_leftparticle_vol_old(0),
    m_rightparticle_vol_old(0),
    m_rightparticle_numPAH(0),
    m_leftparticle_numPAH(0),
    m_children_surf(0),
    m_children_coalescence(0),
    m_Rg(0),
    m_fdim(0),
    m_sqrtLW(0),
    m_LdivW(0),
    m_avg_coalesc(0),
    m_leftchild(NULL),
    m_rightchild(NULL),
    m_parent(NULL),
    m_leftparticle(NULL),
    m_rightparticle(NULL)
{
}


// Initialising constructor.
PAHPrimary::PAHPrimary(real time, const Sweep::ParticleModel &model)
: Primary(time, model),
    m_numcarbon(0),
    m_PAHmass(0),
    m_PAHCollDiameter(0),
    m_numPAH(0),
    m_numprimary(0),
    m_primarydiam(0.0),
    m_children_radius(0),
    m_children_vol(0),
    m_leftparticle_vol_old(0),
    m_rightparticle_vol_old(0),
    m_rightparticle_numPAH(0),
    m_leftparticle_numPAH(0),
    m_children_surf(0),
    m_children_coalescence(0),
    m_Rg(0),
    m_fdim(0),
    m_sqrtLW(0),
    m_LdivW(0),
    m_avg_coalesc(0),
    m_leftchild(NULL),
    m_rightchild(NULL),
    m_parent(NULL),
    m_leftparticle(NULL),
    m_rightparticle(NULL)
{
    // Other parts of the code check for a non-zero composition
    m_comp[0]=1;

    // Select a PAH life story uniformly at random from the database
    Trajectory const *newtrajectory=model.m_PAHDatabase.GetTrajectory(time);
    const int newID = static_cast<int>(newtrajectory->maxID()*rnd());

    AddPAH(time,newID,newtrajectory);
    
    //Get the number of carbon atoms in PAH with id newID from  the database
    UpdatePAHs(time);

    //Update the other properties
    UpdateCache();
}


// Initialising constructor.
PAHPrimary::PAHPrimary(real time, const Sweep::ParticleModel &model, bool noPAH)
: Primary(time, model),
    m_numcarbon(0),
    m_PAHmass(0),
    m_PAHCollDiameter(0),
    m_numPAH(0),
    m_numprimary(0),
    m_primarydiam(0.0),
    m_children_radius(0),
    m_children_vol(0),
    m_leftparticle_vol_old(0),
    m_rightparticle_vol_old(0),
    m_rightparticle_numPAH(0),
    m_leftparticle_numPAH(0),
    m_children_surf(0),
    m_children_coalescence(0),
    m_Rg(0),
    m_fdim(0),
    m_sqrtLW(0),
    m_LdivW(0),
    m_avg_coalesc(0),
    m_leftchild(NULL),
    m_rightchild(NULL),
    m_parent(NULL),
    m_leftparticle(NULL),
    m_rightparticle(NULL)
{
    m_comp[0]=1;

}

/*!
 * Add a PAH to the primary particle
 *
 * @param[in]   time        create time of the PAH
 * @param[in]   ID          integer identifying the PAH in the trajectory database
 * @param[in]   trajectory  Pointer to the trajectory database
*/
void PAHPrimary::AddPAH(real time, int ID, const Trajectory *trajectory)
{
	PAH *new_PAH=new PAH;
	new_PAH->ID=ID;
    new_PAH->m_trajectory=trajectory;
	new_PAH->m_numcarbon=0;
	new_PAH->time_created=time;
	new_PAH->freezetime=0;
	new_PAH->lastupdated=time;
    new_PAH->lastposPAHupdate=0;
	m_PAH.insert(m_PAH.begin(),*new_PAH);
}

// Copy constructor.
PAHPrimary::PAHPrimary(const PAHPrimary &copy) 
{   
    if (this!=&copy)
    {
        *this = copy;
	    if (copy.m_leftchild!=NULL)
	    {
		    CopyTree(&copy);
	    }
    }
}


// Default destructor.
PAHPrimary::~PAHPrimary()
{  
	delete m_leftchild;
	delete m_rightchild;
    // it is not necessary to delete m_leftparticle because 
    // this is also m_leftchild somewhere down the tree  

    // delete the PAH list
    releaseMem();

}

/*!
 * @param[in] source Pointer to the primary to be copied
 * @param[in,out] target Pointer to the target primary
*/
void PAHPrimary::CopyTree( const PAHPrimary *source)
{
    //create the new left and right children 
	PAHPrimary *newleft = new PAHPrimary(source->CreateTime(),*m_pmodel,false);
	PAHPrimary *newright = new PAHPrimary(source->CreateTime(),*m_pmodel,false);
    // copy the properties
	newleft->CopyParts(source->m_leftchild);
	newright->CopyParts(source->m_rightchild);
    //set the pointers to the children and the parent
	m_leftchild=newleft;
	m_rightchild=newright;
    newleft->m_parent=this;
	newright->m_parent=this;

    // the left and right particle are set further down in UpdateAllPointers
	m_leftparticle=NULL;
	m_rightparticle=NULL;
    //copy the subtree
	if (source->m_leftchild->m_leftchild!=NULL) m_leftchild->CopyTree(source->m_leftchild);
	if (source->m_rightchild->m_leftchild!=NULL) m_rightchild->CopyTree(source->m_rightchild);
    //set the leftparticle and rightparticle
	UpdateAllPointers(source);
}
	
PAHPrimary &PAHPrimary::operator=(const Primary &rhs)
{
    if (this != &rhs) {
		const AggModels::PAHPrimary *pahprimary = NULL;
		pahprimary = dynamic_cast<const AggModels::PAHPrimary*>(&rhs);
	//	UpdateAllPointers(pahprimary,this);
    }
    return *this;
}


// Stream-reading constructor.
PAHPrimary::PAHPrimary(std::istream &in, const Sweep::ParticleModel &model)
{
    Deserialize(in, model);
}


//returns a copy of the PAHprimary, the pointers are not corrected. 
// The children are not copied
//It can be necessary to adjust the pointers to the particles afterwards
/*!
 * @param[in] source Pointer to the primary to be copied
 * @param[in,out] target Pointer to the target primary
*/
void PAHPrimary::CopyParts( const PAHPrimary *source)
{
	SetCollDiameter(source->CollDiameter());
	SetMobDiameter(source->MobDiameter());
	SetSphDiameter(source->SphDiameter());
	m_PAHCollDiameter=source->m_PAHCollDiameter;
	SetSurfaceArea(source->SurfaceArea());
	m_time=source->m_time;
	m_PAHmass=source->m_PAHmass;
    m_leftchild=source->m_leftchild;               
	m_rightchild=source->m_rightchild;
	m_leftparticle=source->m_leftparticle;
	m_rightparticle=source->m_rightparticle;
	m_parent=source->m_parent;
	SetMass(source->Mass());
	m_numPAH=source->m_numPAH;
	m_primarydiam=source->m_primarydiam;
	m_sqrtLW=source->m_sqrtLW;
	m_LdivW=source->m_LdivW;
    m_pmodel=source->m_pmodel;
    m_surf=source->m_surf;
    m_vol=source->m_vol;
    m_children_surf=source->m_children_surf;
    m_children_vol=source->m_children_vol;
    m_children_radius=source->m_children_radius;
    m_children_coalescence=source->m_children_coalescence;
    m_rightparticle_numPAH=source->m_rightparticle_numPAH;
    m_leftparticle_numPAH=source->m_leftparticle_numPAH;
    m_leftparticle_vol_old=source->m_leftparticle_vol_old;
    m_rightparticle_vol_old=source->m_rightparticle_vol_old;
    m_fdim=source->m_fdim;
    m_Rg=source->m_Rg;
    m_avg_coalesc=source->m_avg_coalesc;
    m_numcarbon=source->m_numcarbon;
	vector<PAH>::const_iterator j;
	//for (j=source->m_PAH.begin(); j!=source->m_PAH.end(); ++j) {
	//	m_PAH.insert(m_PAH.end(),PAH(*j));
	//}
    m_PAH.clear();
    m_PAH.insert(m_PAH.end(),source->m_PAH.begin(),source->m_PAH.end());
	
}

//selects a uniformly chosen primary  
/*!
 * @param[in] root Pointer to the primary root node
*/
PAHPrimary *PAHPrimary::SelectRandomSubparticle()
{
	double ran=rnd();
    int target=int(ran*(m_numprimary));
	return SelectRandomSubparticleLoop(target);

}
/*!
 * @param[in] source Pointer to the primary root node to be copied
 * @param[in] target The primary to be selected
*/
PAHPrimary *PAHPrimary::SelectRandomSubparticleLoop(int target)
{	
	if (m_leftchild==NULL) return this;
	if (target<=m_leftchild->m_numprimary)
	{
		return m_leftchild->SelectRandomSubparticleLoop(target);
	}
	else
	{
		return m_rightchild->SelectRandomSubparticleLoop(target-(m_leftchild->m_numprimary));
	}
}


// This function is used to find the the left and right particles when the entire particle tree is duplicated
/*!
 * @param[in] source Pointer to the primary to be copied
 * @param[in,out] target Pointer to the target primary
*/
void PAHPrimary::UpdateAllPointers( const PAHPrimary *original)
{	
    //the primary has no children => there are no left and right particles
    if (original->m_leftchild == NULL) {
		m_leftparticle=NULL;
		m_rightparticle=NULL;			
    } else {	
		int depth=0;
		int i;
		bool path[1000];
		PAHPrimary *orgparticle;
        //update the left particle 
		orgparticle=original->m_leftparticle;
        // search for the left particle in the original tree and store the pathway in path
		while (orgparticle!=original)
		{
			if (orgparticle->m_parent->m_leftchild==orgparticle)
				path[depth]=true;
			else
				path[depth]=false;
			orgparticle=orgparticle->m_parent;
			depth++;
		}
        // go the same path in the target tree
		m_leftparticle=m_leftchild;
		for (i=depth-2;i>-1;i--)
		{
			if (path[i]==true)
            {
				m_leftparticle=m_leftparticle->m_leftchild;
            }
			else 
            {
				m_leftparticle=m_leftparticle->m_rightchild;
            }
		}
		
         //update the right particle 
		depth=0;
		orgparticle=original->m_rightparticle;
        // search for the right particle in the original tree and store the pathway in path
		while (orgparticle!=original)
		{
			if (orgparticle->m_parent->m_leftchild==orgparticle)
				path[depth]=true;
			else
				path[depth]=false;
			depth++;
			orgparticle=orgparticle->m_parent;
		}
        // go the same path in the target tree
		m_rightparticle=m_rightchild;
		for (i=depth-2;i>-1;i--)
		{
			if (path[i]==true)
            {
				m_rightparticle=m_rightparticle->m_leftchild;
            }
			else 
            {
				m_rightparticle=m_rightparticle->m_rightchild;
            }
		}
	}
}


// coagulates this particle and rhs
/*!
 * @param[in] rhs Pointer to the particle to be coagulate with this particle
*/
PAHPrimary &PAHPrimary::Coagulate(const Primary &rhs)
{
    vector<PAH>::const_iterator j;
	const PAHPrimary *rhsparticle = NULL;
	rhsparticle = dynamic_cast<const AggModels::PAHPrimary*>(&rhs);

	//only one PAH in rhs or this particle -> condensation or inception process. 
	if ( (rhsparticle->m_numPAH==1 || this->m_numPAH==1 ))
	{
        if (rhsparticle->Numprimary()>1)
        {
            PAHPrimary copy_rhs(*rhsparticle);
            PAHPrimary *target = copy_rhs.SelectRandomSubparticle();
			//for (j=m_PAH.begin(); j!=m_PAH.end(); ++j) {
			//	target->m_PAH.insert(target->m_PAH.end(),PAH(*j));
            //}
            target->m_PAH.insert(target->m_PAH.end(),m_PAH.begin(),m_PAH.end());
            target->UpdatePrimary();
            CopyParts(&copy_rhs);
            CopyTree(&copy_rhs);
            copy_rhs.m_leftchild=NULL;
	        copy_rhs.m_rightchild=NULL;
	        copy_rhs.m_parent=NULL;
            copy_rhs.m_leftparticle=NULL;
            copy_rhs.m_rightparticle=NULL;
        }
        else
        {
            // particle has more then one primary select the primary where 
            // the PAH condenses to and add it to the list
		    if (m_leftchild!=NULL)
		    {
			    PAHPrimary *target = SelectRandomSubparticle();
	    //		for (j=rhsparticle->m_PAH.begin(); j!=rhsparticle->m_PAH.end(); ++j) {
	    //			target->m_PAH.insert(target->m_PAH.end(),PAH(*j));
	    //		}
                target->m_PAH.insert(target->m_PAH.end(),rhsparticle->m_PAH.begin(),rhsparticle->m_PAH.end());
                target->UpdatePrimary();

		    }
		    else
		    {
		    //	for (j=rhsparticle->m_PAH.begin(); j!=rhsparticle->m_PAH.end(); ++j) {
		    //		m_PAH.insert(m_PAH.end(),PAH(*j));
		    //	}
                m_PAH.insert(m_PAH.end(),rhsparticle->m_PAH.begin(),rhsparticle->m_PAH.end());
                UpdatePrimary();
		    }
        }
		UpdateCache();
        //Check the coalescence ratio 
        CheckCoalescence();

	}
	
    else
    {   
            //coagulation process 
            PAHPrimary *newleft = new PAHPrimary;
		    PAHPrimary *newright = new PAHPrimary;
            PAHPrimary copy_rhs(*rhsparticle);
		    //rhsparticle = dynamic_cast<const AggModels::PAHPrimary*>(&rhs);

           // bool print=false;
           // if (this->m_numprimary>2 && notconstpah->m_numprimary>2)
           // {
            //    PrintTree("before1");
           //     notconstpah->PrintTree("before2");
             //   print=true;
              //  cout << "printing tree"<<endl;
            //}

            //select where to add the second particle
		    if (rnd()>0.5)
		    {
			    newleft->CopyParts(this);
			    newright->CopyParts(&copy_rhs);
		    }
		    else
		    {
			    newright->CopyParts(this);
			    newleft->CopyParts(&copy_rhs);
		    }
            //set the pointers
		    m_leftchild=newleft;
		    m_rightchild=newright;
		    newright->m_parent=this;
		    newleft->m_parent=this;
            //set the pointers to the parent node
		    if (newleft->m_leftchild!=NULL)
		    {
			    newleft->m_leftchild->m_parent=newleft;
			    newleft->m_rightchild->m_parent=newleft;
		    }
		    if (newright->m_leftchild!=NULL)
		    {
			    newright->m_leftchild->m_parent=newright;
			    newright->m_rightchild->m_parent=newright;
		    }
            m_children_coalescence=0;
		    UpdateCache();
            //select the primaries that are touching
		    this->m_leftparticle=m_leftchild->SelectRandomSubparticle();
		    this->m_rightparticle=m_rightchild->SelectRandomSubparticle();

            //initialise the variables used to calculate the coalesence ratio
            m_children_vol=m_leftparticle->m_vol+m_rightparticle->m_vol;
            m_children_surf=(m_leftparticle->m_surf+m_rightparticle->m_surf);
            m_leftparticle_vol_old=m_leftparticle->m_vol;
            m_rightparticle_vol_old=m_rightparticle->m_vol;
            m_leftparticle_numPAH=m_leftparticle->m_numPAH;
            m_rightparticle_numPAH=m_rightparticle->m_numPAH;
            m_children_radius=pow(3.0/(4.0*PI)*(m_children_vol),(ONE_THIRD));
            m_children_coalescence=CoalescenceLevel();  
            CheckCoalescence();
            //must set all the pointer to NULL otherwise the delete function 
            //will also delete the children 
	        copy_rhs.m_leftchild=NULL;
	        copy_rhs.m_rightchild=NULL;
	        copy_rhs.m_parent=NULL;
            copy_rhs.m_leftparticle=NULL;
            copy_rhs.m_rightparticle=NULL;

        //    if (fabs(surfacebeforerhs+surfacebefore-m_surf)/m_surf > 1e-6)
        //    {
        //        cout << "error" << surfacebeforerhs<<' ' <<surfacebefore << ' '<<m_surf<<endl;
        ////         PrintTree("after");
        //    }
        // if (print)
    //     PrintTree("after");
        
	}
    return *this;
}

//returns the CoalescenceLevel of the two primaries that are connected by this node
double PAHPrimary::CoalescenceLevel() 
{
    if (m_leftparticle!=NULL)
    {

        double dV=0;
        
        
        //make sure that the volume growth does not arise from a former coalescence event 
        //(number of PAHs must remain const, or increase by 1 in case of a condensation
        // event
        // if the volume changes due to a coalescence event the chilren surface is changed in 
        // PAHPrimary::ChangePointer(PAHPrimary *source, PAHPrimary *target)
        
        // m_leftparticle->m_numPAH is the actual number of PAHs in the left one of the two
        // touching primary particles.

        // m_leftparticle_numPAH is the the number of PAHs in the left one of the two touching
        // particles last time this method was called.

        // The following test checks that at most one PAH has been added to the left one
        // of the two touching primaries since this method was last called.
        if (m_leftparticle->m_numPAH - m_leftparticle_numPAH<2) 
        {   
            //calculate the volume change of the left particle in the last timestep
            dV+= m_leftparticle->m_vol - m_leftparticle_vol_old;         
        }

        //adjust the number of PAHs and vol for the next step
        m_leftparticle_numPAH=m_leftparticle->m_numPAH;
        m_leftparticle_vol_old=m_leftparticle->m_vol;

        //make sure that the volume growth does not arise from a coalescence event 
        //(number of PAHs must remain const, or increase by1 in case of a condensation
        // event
        if (m_rightparticle->m_numPAH-m_rightparticle_numPAH<2)   
        {   
            //calculate the volume change of the right particle in the last timestep
            dV+= m_rightparticle->m_vol - m_rightparticle_vol_old;          
        }
        
        //adjust the number of PAHs and vol for the next step
        m_rightparticle_numPAH=m_rightparticle->m_numPAH;
        m_rightparticle_vol_old=m_rightparticle->m_vol;

        //update the children volume, this value is used to calculate dV n the next timestep
        m_children_vol=m_rightparticle->m_vol+m_leftparticle->m_vol;
        //calculate dS, at the moment it is assumed that the particles always grow
        double ct=m_pmodel->Components(0)->CoalescThresh();
        const double dS=dV*ct/m_children_radius;
        //update the radius for the next event
        m_children_radius=pow(3.0/(4.0*PI)*(m_children_vol),(ONE_THIRD));
        m_children_surf+=dS;

        const double spherical_surface=4*PI*m_children_radius*m_children_radius;
       // double two_1_3=pow(2,-1*ONE_THIRD);
        const double two_1_3=0.79370052231642452;
        double clevel= ((spherical_surface/m_children_surf)-two_1_3)/(1-two_1_3);
        if (clevel<0)
            return 0;
        else
            return clevel;
    }
    else 
        return 0;
}

//calculates the fractal dimension of the particle and stores it in m_fdim
void PAHPrimary::CalcFractalDimension()
{
	Sweep::Imaging::ParticleImage img;
    // construct the particle by colliding the primary particles
	img.constructSubParttree(this);
	double L,W;
    // calculate the length and the width of the particle 
    img.LengthWidth(L,W);
    // calculate the radius of gyration
    m_Rg=img.RadiusofGyration();
	m_sqrtLW=sqrt(L*W);
	m_LdivW=L/W;
    m_Rg=m_Rg*1e-9;
    m_fdim=log((double)m_numprimary)/log(2*m_Rg/(m_primarydiam/m_numprimary));
  /*  if (m_fdim>0 && m_fdim<3)
    {    
       string filename;
       filename=cstr(m_fdim)+".3d";
       ofstream out;
       out.open(filename.c_str());
       img.Write3dout(out,0,0,0);
       out.close();
    }*/
    
}   

// sets all the childrenproperties to zero, this function is used after the children are coalesced 
void PAHPrimary::ResetChildrenProperties()
{
            m_children_coalescence=0.0;
            m_children_surf=0.0;
            m_children_vol=0.0;
            m_rightparticle_vol_old=0.0;
            m_leftparticle_vol_old=0.0;
            m_leftparticle_numPAH=0;
            m_rightparticle_numPAH=0;
            m_avg_coalesc=0;
}

//merges the left and the right primary particle to one primary particle
PAHPrimary &PAHPrimary::Merge()
{
//    if(this->m_numprimary>5)
 //       cout<<"tset";
 //      PrintTree("before.inp");

    vector<PAH>::const_iterator j;
      //make sure this primary has children to merge
	  if(m_leftchild!=NULL)
           {
		if ( m_leftchild==m_leftparticle && m_rightchild==m_rightparticle)
		{	
            //this node has only two primaries in its subtree
            //it is possible that this node is not the root node and belongs to a bigger particle
            //copy the PAHs of both children to the parent node

            //m_PAH is empty, therefore no need to append 
            m_PAH=m_rightparticle->m_PAH;

            m_PAH.insert(this->m_PAH.end(),m_leftparticle->m_PAH.begin(),m_leftparticle->m_PAH.end());
            //update the pointers that pointed to the two former children
			ChangePointer(m_leftchild,this);
			ChangePointer(m_rightchild,this);
            //delete the children
			m_leftchild->ReleaseMem();
			m_rightchild->ReleaseMem();
            m_leftchild->m_parent=NULL;
            m_rightchild->m_parent=NULL;
            delete m_leftchild;
            delete m_rightchild;                    
			m_leftchild=NULL;
			m_rightchild=NULL;
			m_leftparticle=NULL;
			m_rightparticle=NULL;
            //set the children properties to zero, this node has no more children
            ResetChildrenProperties();
            UpdatePrimary();
		}


		else
		{
			if (m_leftchild->m_numprimary<m_rightchild->m_numprimary)
			{
				//append to left subtree because there are fewer primaries 
                //this is only to keep the tree balanced
				PAHPrimary *oldleftparticle=m_leftparticle;

                //copy the PAHs

				//for (j=oldleftparticle->m_PAH.begin(); j!=oldleftparticle->m_PAH.end(); ++j) {
				//	m_rightparticle->m_PAH.insert(m_rightparticle->m_PAH.end(),PAH(*j));
				//}
                m_rightparticle->m_PAH.insert(m_rightparticle->m_PAH.end(),oldleftparticle->m_PAH.begin(),oldleftparticle->m_PAH.end());
                m_rightparticle->UpdatePrimary();
                //set the pointers from the leftprimary to the rightprimary
                //this will be the new bigger primary
				oldleftparticle->ChangePointer(oldleftparticle,m_rightparticle);
                m_rightparticle->ChangePointer(m_rightparticle,m_rightparticle);
                //set the pointer to the parent node 
				if (oldleftparticle->m_parent->m_leftchild==oldleftparticle)
				{
					oldleftparticle->m_parent->m_leftchild=m_rightchild;
				}
				else
				{
					oldleftparticle->m_parent->m_rightchild=m_rightchild;
				}
				m_rightchild->m_parent=oldleftparticle->m_parent;

			
				PAHPrimary *oldleftchild=m_leftchild;
				PAHPrimary *oldparent=m_parent;
                //copy the properties of the former leftchild to this node
				CopyParts(oldleftchild);
				m_parent=oldparent;
                if (m_leftchild!=NULL)
                {
				    m_rightchild->m_parent=this;
				    m_leftchild->m_parent=this; 
                } 
				oldleftchild->ReleaseMem();

			}

			else
			{
				//append to right subtree
				PAHPrimary *oldrightparticle=m_rightparticle;
				
			//	for (j=oldrightparticle->m_PAH.begin(); j!=oldrightparticle->m_PAH.end(); ++j) {
			//		m_leftparticle->m_PAH.insert(m_leftparticle->m_PAH.end(),PAH(*j));
			//	}
                m_leftparticle->m_PAH.insert(m_leftparticle->m_PAH.end(),oldrightparticle->m_PAH.begin(),oldrightparticle->m_PAH.end());
                m_leftparticle->UpdatePrimary();
                oldrightparticle->ChangePointer(oldrightparticle,m_leftparticle);
                m_leftparticle->ChangePointer(m_leftparticle,m_leftparticle);
				if (oldrightparticle->m_parent->m_leftchild==oldrightparticle)
				{
					oldrightparticle->m_parent->m_leftchild=m_leftchild;
				}
				else
				{
					oldrightparticle->m_parent->m_rightchild=m_leftchild;
				}
				m_leftchild->m_parent=oldrightparticle->m_parent;

				//ReleaseMem(oldrightparticle);
				PAHPrimary *oldrightchild=m_rightchild;
				PAHPrimary *oldparent=m_parent;
				CopyParts(oldrightchild);
				m_parent=oldparent;
                if (m_leftchild!=NULL)
                {
				    m_rightchild->m_parent=this;
				    m_leftchild->m_parent=this; 
                } 
				oldrightchild->ReleaseMem();
                
			}

		}
         
        UpdateCache();
  //      PrintTree("after.inp");
        }
		return *this;
}

void PAHPrimary::ReleaseMem()
{
    m_PAH.clear();  
}
/*!//add explanation
 * @param[in] source Pointer to the original particle 
 * @param[in,out] target Pointer to the new particle
*/
void PAHPrimary::ChangePointer(PAHPrimary *source, PAHPrimary *target)
{    
		if(m_rightparticle==source){
			m_rightparticle=target;
            double sphericalsurface=
                4*PI*pow(3*(m_leftparticle->Volume()+m_rightparticle->Volume())/(4*PI),TWO_THIRDS);
            m_children_surf=sphericalsurface/(m_children_coalescence*0.2063+0.7937);    //sphericalsurface/(m_children_coalescence*(1-2^(-1/3))+2^(-1/3))
		}
		if(m_leftparticle==source){
			m_leftparticle=target;
            double sphericalsurface=
                4*PI*pow(3*(m_leftparticle->Volume()+m_rightparticle->Volume())/(4*PI),TWO_THIRDS);
            m_children_surf=sphericalsurface/(m_children_coalescence*0.2063+0.7937);    //sphericalsurface/(m_children_coalescence*(1-2^(-1/3))+2^(-1/3))

		}

    // Update the tree above this sub-particle.
    if (m_parent != NULL) {
        m_parent->ChangePointer(source,target);
    }

}

/*!
 * @param[in] t Time
*/
void PAHPrimary::UpdatePAHs(double t)
{
	if (m_leftchild!=NULL)
	{
		m_leftchild->UpdatePAHs(t);
		m_rightchild->UpdatePAHs(t);
	}
    else
    {   
        unsigned int j=0;
	    for (vector<PAH>::iterator i=m_PAH.begin(); i!=m_PAH.end(); ++i) {
		    //for (j=0;j<i->time.size();j++)
            double minPAH=m_pmodel->Components(0)->MinPAH();
            if (m_numPAH>=minPAH)
            {
                double growthfact=m_pmodel->Components(0)->GrowthFact();
                double deltat=t-i->lastupdated;
                i->freezetime=i->freezetime+deltat*(1.0-growthfact);
            }
            double seektime=i->m_trajectory->StartTime()+t-(i->time_created)-(i->freezetime);
		    for (j=i->lastposPAHupdate;j<i->m_trajectory->alltrajectories.at(i->ID).time.size();j++)
            {  
                //if (this->m_numPAH>1000) i->freezetime+=rnd()*(t-i->lastupdated);
			    //i->m_numcarbon=i->n_carbon_t.at(j);
			    i->m_numcarbon=i->m_trajectory->alltrajectories.at(i->ID).n_carbon_t.at(j);
			    //if (i->time.at(j)>=t-i->time_created-i->freezetime)
                
			    if (i->m_trajectory->alltrajectories.at(i->ID).time.at(j)>=0.999*seektime)
			    {
                    i->lastposPAHupdate=j;
				    break;            
			    }
		    }
		    i->lastupdated=t;
        }

        UpdatePrimary();
	}
}


void PAHPrimary::UpdateCache(void)
{
    UpdateCache(this);
}


bool PAHPrimary::CheckCoalescence()
{
    bool hascoalesced=false;
   // cout << m_children_coalescence<<endl;
 //  if (m_children_coalescence> m_pmodel->Components(0)->CoalescThresh() && m_leftparticle!=NULL)
    if (m_children_coalescence> 0.99 && m_leftparticle!=NULL)
        {
           // PrintTree("before.inp");
          //   cout <<"merging"<<m_children_coalescence<<endl;
             Merge();
           //  PrintTree("after.inp");
           hascoalesced=true;
           //check again because this node has changed
         CheckCoalescence();
        }
    if (m_leftchild!=NULL)
    {
        hascoalesced=m_leftchild->CheckCoalescence();
        hascoalesced=m_rightchild->CheckCoalescence();
    }
    UpdateCache();
    return hascoalesced;
}


//this function updates a primary particle
void PAHPrimary::UpdatePrimary(void)
{	 
    double sootdens=0;
	m_numcarbon=0;
    m_PAHmass=0;
	m_PAHCollDiameter=0;
	m_numPAH=0;
    unsigned int maxcarbon=0;
    for (vector<PAH>::iterator i=m_PAH.begin(); i!=m_PAH.end(); ++i) {
        m_numcarbon += i->m_numcarbon;
		maxcarbon=max(maxcarbon,i->m_numcarbon);    // search for the largest PAH in the PRimary, in Angstrom
		m_numPAH++;
    }
	m_PAHmass=m_numcarbon*1.9945e-26;        //convert to kg, hydrogen atoms are not considered
    m_PAHCollDiameter=sqrt(maxcarbon*2.0/3.);
	m_PAHCollDiameter*=2.4162*1e-10;				 //convert from Angstrom to m

    if(m_pmodel->ComponentCount()!=1)        //at the moment we have only one component: soot
    {
        cout<<"Warning: Model contains more then one component. At the moment only soot is supported"<<endl;
    }
    else
    {   
       sootdens=m_pmodel->Components(0)->Density();        //in kg/m^3
    }
	 m_mass=m_PAHmass;
	 m_vol=m_PAHmass/sootdens;
	 m_diam = pow(6.0 * m_vol / PI, ONE_THIRD);    
     m_dmob = m_diam;
     m_dcol = max(m_diam,m_PAHCollDiameter);
     m_surf = PI * m_diam * m_diam;
     m_primarydiam = m_diam;
     m_avg_coalesc=0;
}



void PAHPrimary::Reset()
{
	m_numcarbon=0;
	m_primarydiam=0.0;
    m_surf=0;
    m_vol=0;
    m_PAH.clear();
    m_avg_coalesc=0;
}
/*!
 * @param[in] root The root node of this particle
*/
void PAHPrimary::UpdateCache(PAHPrimary *root)
{	
    //Update the children 
	if (m_leftchild!=NULL)
	{
		m_leftchild->UpdateCache(root);
		m_rightchild->UpdateCache(root);
		m_numprimary=m_leftchild->m_numprimary+m_rightchild->m_numprimary;
	}
    //this is a primary and the number of primaries below this node is one (this node and no children)
	else
	{
		m_numprimary=1;
        m_avg_coalesc=0;
	}
     
    //this is not a primary, sum up the properties 
	if (m_leftchild!=NULL) 
    {
        // remove the PAHs from this node to free memory
        Reset();
		m_surf = m_leftchild->m_surf+m_rightchild->m_surf;
		m_primarydiam = (m_leftchild->m_primarydiam+m_rightchild->m_primarydiam);
        m_vol=m_leftchild->m_vol+m_rightchild->m_vol;
        m_numPAH = m_leftchild->m_numPAH+m_rightchild->m_numPAH;
        m_primarydiam = (m_leftchild->m_primarydiam+m_rightchild->m_primarydiam);
        m_mass=(m_leftchild->m_mass+m_rightchild->m_mass);
        m_PAHCollDiameter=max(m_leftchild->m_PAHCollDiameter,m_rightchild->m_PAHCollDiameter);
        m_numcarbon=m_leftchild->m_numcarbon+m_rightchild->m_numcarbon;
        // calculate the coalescence level of the two primaries connected by this node
        m_children_coalescence=CoalescenceLevel();
        if (m_children_coalescence>1)
            m_children_coalescence=1;
        //sum up the avg coal level 
        m_avg_coalesc=m_children_coalescence+m_leftchild->m_avg_coalesc+m_rightchild->m_avg_coalesc;
        
        // calculate the different diameters only for the root node because this goes into the 
        // particle tree and gets used by the coagulation kernel
        if (this==root)
        {
             //spherical eqiv radius
            double spherical_radius=pow(3*m_vol/(4*PI),ONE_THIRD);         
            m_diam=2*spherical_radius;
            // there are m_numprimary-1 connections between the primary particles          
            m_avg_coalesc=m_avg_coalesc/(m_numprimary-1);
            //approxmiate the surface of the particle
            const real numprim_1_3=pow(m_numprimary,-0.333333);
            m_surf=4*PI*spherical_radius*spherical_radius/
                (m_avg_coalesc*(1-numprim_1_3)+numprim_1_3);

            //calculate the surface equivalent radius
           // const double radius_surf=sqrt(m_surf/(4*PI));
            // the average between the surface and voluem equiv diameter
            //const double meandiam=spherical_radius+radius_surf;            //2*0.5*(radius(vol)+radius(sphere))
            const double aggcolldiam=(6*m_vol/m_surf)*
                                     pow(pow(m_surf,3)/(36*PI*m_vol*m_vol),(1.0/1.8));
            // the maximum of the largest PAH diameter and  
            // the average between the surface and voluem equiv diameter       
            const double cdiam=max(aggcolldiam,m_PAHCollDiameter);      
            m_dmob = aggcolldiam;
            SetCollDiameter(cdiam);
        }
        else 
        {
            m_diam=0;
            m_dmob=0;
        }
    }
}

/*!
 * @param[in] filename Output filename
*/
void PAHPrimary::PrintTree(string filename)
{
  ofstream out;
  out.open(filename.c_str());
  out << "digraph unix {"<<endl;
  out <<"graph [rankdir = \"LR\"];"<<endl;
  PrintTreeLoop(out);
  out << "}"<<endl;
  out.close();
}

/*!
 * @param[in] out Output stream
*/
void PAHPrimary::PrintTreeLoop(std::ostream &out)
{ 
  if (m_leftchild!=NULL)
  { //out<<"leftchild "<<10E8*m_leftchild->SphDiameter()<<endl;   
	//out<<"rightchild "<<10E8*m_rightchild->SphDiameter()<<endl;    
      out<<"\" "<<this<<"\" "<<" [shape = \"record\" label = \"surf="<<this->m_surf<<"|m_children_surf="<<this->m_children_surf<<"|m_vol="<<this->m_vol<<"|"<<this->m_children_coalescence<<"|"<<this<<"\"];"<<endl;
	out<<"\" "<<this->m_leftchild<<"\" "<<" [shape = \"record\" label = \"surf="<<this->m_surf<<"|m_children_surf="<<this->m_children_surf<<"|m_vol="<<this->m_vol<<"|"<<this->m_children_coalescence<<"|"<<this<<"\"];"<<endl;
	out<<"\" "<<this->m_rightchild<<"\" "<<" [shape = \"record\" label = \"surf="<<this->m_surf<<"|m_children_surf="<<this->m_children_surf<<"|m_vol="<<this->m_vol<<"|"<<this->m_children_coalescence<<"|"<<this<<"\"];"<<endl;
	out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_leftchild<<"\"; "<<endl;
	out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_rightchild<<"\"; "<<endl;
	out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_leftparticle<<"\"[label=\""<<this<<"\",color=\"blue\"]; "<<endl;
	out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_rightparticle<<"\"[label=\""<<this<<"\",color=\"blue\"]; "<<endl;
	m_leftchild->PrintTreeLoop(out);
    m_rightchild->PrintTreeLoop(out);
  }

  else
  {  
      //out<<"\" "<<this<<"\" "<<" [label = \""<<2*pow((this->vol_sinter)*3/(4*PI),ONE_THIRD)<<"\"];"<<endl;
      out<<"\" "<<this<<"\" "<<" [shape = \"record\" label = \"surf="<<this->m_surf<<"|m_children_surf="<<this->m_children_surf<<"|m_vol="<<this->m_vol<<"|"<<this->m_children_coalescence<<"|"<<this<<"\"];"<<endl;
  }
}


// Creates an aggregation data cache for this primary type.
AggModels::PAHCache *const PAHPrimary::CreateAggCache(ParticleCache &pcache) const
{
    PAHCache *cache = 
        static_cast<PAHCache*>(ModelFactory::CreateAggCache(AggModels::PAH_ID, pcache));
    if (cache != NULL) *cache = *this;
    return cache;
}


 const PAHPrimary *PAHPrimary::RightChild() const
{
    return m_rightchild;
}

 const PAHPrimary *PAHPrimary::LeftChild() const
{
    return m_leftchild;
}

double PAHPrimary::PAHCollDiameter() const
{
    return m_PAHCollDiameter;
}

double PAHPrimary::Rg() const
{
    return m_Rg;
}

double PAHPrimary::Fdim() const
{
    return m_fdim;
}

double PAHPrimary::PrimaryDiam() const
{
    return m_primarydiam;
}

double PAHPrimary::LdivW() const
{
    return m_LdivW;
}

int PAHPrimary::Numprimary() const
{
    return m_numprimary;
}

int PAHPrimary::NumCarbon() const
{
    return m_numcarbon;
}

int PAHPrimary::NumPAH() const
{
    return m_numPAH;
}

double PAHPrimary::sqrtLW() const
{
    return m_sqrtLW;
}

double PAHPrimary::AvgCoalesc() const
{
    return m_avg_coalesc;
}

// READ/WRITE/COPY.

// Returns a copy of the model data.
PAHPrimary *const PAHPrimary::Clone(void) const
{
    return new PAHPrimary(*this);
}


// AGGREGATION MODEL.

// Returns the aggregation model which this primary describes.
AggModels::AggModelType PAHPrimary::AggID(void) const {return AggModels::PAH_ID;}



// Writes the object to a binary stream.
void PAHPrimary::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

		double val = 0.0;
		// Write number of PAHs.
        val = (double) m_numPAH;
        out.write((char*)&val, sizeof(val));

		// Write PAHcollisiondiameter
        val = (double) m_PAHCollDiameter;
        out.write((char*)&val, sizeof(val));

	    val = (double)m_numprimary;
        out.write((char*)&val, sizeof(val));

	    val = (double)m_sqrtLW;
        out.write((char*)&val, sizeof(val));

	    val = (double)m_LdivW;
        out.write((char*)&val, sizeof(val));

	    val = (double)m_primarydiam;
        out.write((char*)&val, sizeof(val));

        val = (double)m_fdim;
        out.write((char*)&val, sizeof(val));

        val = (double)m_Rg;
        out.write((char*)&val, sizeof(val));

		// write the PAH stack
		/*PAH currPAH;
		for (int i=0; i!=m_numPAH; ++i) {
			currPAH = m_PAH[i];
                        out.write((char*)&currPAH.ID, sizeof(currPAH.ID));
                        out.write((char*)&currPAH.time_created, sizeof(currPAH.time_created));
                        out.write((char*)&currPAH.m_numcarbon, sizeof(currPAH.m_numcarbon));
                }
*/
		// Write PAHmass
        val = (double) m_PAHmass;
        out.write((char*)&val, sizeof(val));
		


        // Output base class.
        Primary::Serialize(out);

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, PAHPrimary::Serialize).");
    }
}

// Reads the object from a binary stream.
void PAHPrimary::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

		double val = 0.0;
		// Read number of PAHs.
        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_numPAH = (int)val;
        
		// Read PAHcolldiamter.
        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_PAHCollDiameter = (real)val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_numprimary = (int)val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_sqrtLW = (real)val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_LdivW = (real)val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_primarydiam = (real)val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_fdim = (real)val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_Rg = (real)val;	


	    // Read PAHs.
	/*	
		for (int i=0; i!=m_numPAH; ++i) {
			PAH currPAH;
			//currPAH = new PAH;
			in.read(reinterpret_cast<char*>(&currPAH.ID), sizeof(currPAH.ID));
			in.read(reinterpret_cast<char*>(&currPAH.time_created), sizeof(currPAH.time_created));
			in.read(reinterpret_cast<char*>(&currPAH.m_numcarbon), sizeof(currPAH.m_numcarbon));
			m_PAH.push_back(currPAH); 
		}
*/
		// Read PAHmass.
        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_PAHmass = (real)val;
		m_leftchild=NULL;
		m_rightchild=NULL;
		m_parent=NULL;
		m_leftparticle=NULL;
		m_rightparticle=NULL;
        switch (version) {
            case 0:
                // Read base class.
                Primary::Deserialize(in, model);

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, PAHPrimary::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, PAHPrimary::Deserialize).");
    }
}
