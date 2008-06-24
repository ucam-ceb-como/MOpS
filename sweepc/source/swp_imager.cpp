#include "swp_imager.h"
#include "rng.h"
#include "string_functions.h"
#include <fstream>
#include <vector>

using namespace Sweep;
using namespace Sweep::Imaging;
using namespace std;
using namespace Strings;

const real ParticleImage::m_necking = 1.15;

// CONSTRUCTORS AND DESTRUCTORS.

// Default contructor.
ParticleImage::ParticleImage(void)
{
}

// Initialising constructor.
ParticleImage::ParticleImage(const Particle &sp)
{
    Construct(sp);
}

// Destructor.
ParticleImage::~ParticleImage(void)
{
}

// PARTICLE IMAGE DATA CONSTRUCTION.

// Constructs the particle image from the given particle.
void ParticleImage::Construct(const Particle &sp)
{
    switch(m_creg) {
        case FreeMol:
            // Free-molecular is the default regime.
        default:
            constructAgg_FM(sp);
    }
}


// RENDERING FUNCTIONS.

// Draws the particle image to a POVRAY file.
void ParticleImage::WritePOVRAY(std::ofstream &file)
{
    if (file.good()) {
        string line;
        real val = 0.0;

        // vector of arrays to store primary coordinates.  First
        // 3 values are the cartesian coordinates, final value
        // is the primary radius.
        vector<fvector> coords;

        // Write ParticleDiameter argument to POV file.
        val  = m_root.Radius() * 2.0;
        line = "#declare ParticleDiameter = " + cstr(val) + ";\n";
        file.write(line.c_str(), line.length());

        // Write aggregate opening declaration.
        line = "#declare MyParticle = blob {\n";
        file.write(line.c_str(), line.length());

        // Write threshold radius based on necking parameter.
        val  = pow(1.0 - (1.0/(m_necking*m_necking)), 2.0);
        line = "  threshold " + cstr(val) + "\n";
        file.write(line.c_str(), line.length());

        // Get the primary coordinates from the aggregate tree.
        m_root.GetPriCoords(coords);

        // Write the primaries to the POV-RAY file.
        for (unsigned int i=0; i!=coords.size(); ++i) {
            val  = coords[i][3] * m_necking;
            line = "sphere {<" + cstr(coords[i][0]) + ", " + cstr(coords[i][1]) + 
                   ", " + cstr(coords[i][2]) + ">, " + cstr(val) + ", 1.0}\n";
            file.write(line.c_str(), line.length());
        }

        // Write blob translation.

        // Write closing brace for MyParticle declaration.
        line = "}\n";
        file.write(line.c_str(), line.length());

        // Write include command for tem_single.pov, which defines TEM style
        // for a single particle.
        line = "#include \"particle.pov\"\n";
        file.write(line.c_str(), line.length());
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ParticleImage::WritePOVRAY).");
    }
}


// AGGREGATE SPHERE-TREE CONSTRUCTORS (FREE-MOLECULAR).

// Constructs a PNode sphere-tree aggregate from the given 
// particle using free-molecular collision dynamics.
void ParticleImage::constructAgg_FM(const Particle &sp)
{
    // Clear the current image data structure.
    m_root.Clear();

    // Need to determine the type of particle this is.
    if (sp.ParticleModel()->UseSubPartTree()) {
        // This particle uses the sub-particle tree, so 
        // we can usethat to construct the sphere-tree for
        // image output.
        // TODO:  Complete sub-particle tree TEM output.
    } else {
        const AggModels::SurfVolPrimary *svp = NULL;
        const AggModels::PriPartPrimary *ppp = NULL;

        switch(sp.ParticleModel()->AggModel()) {
            case AggModels::Spherical_ID:
                // Spherical particle model, just draw
                // a single sphere.
                m_root.Insert(sp.SphDiameter() * 0.5e9); // Convert to nm.
                return;
            case AggModels::SurfVol_ID:
                // Surface-volume model, construct a tree
                // of identical primaries, estimating the primary
                // count and diameter.
                svp = dynamic_cast<const AggModels::SurfVolPrimary*>(sp.Primary());
                if (svp != NULL) uniformAgg_FM(svp->PP_Count(), svp->PP_Diameter());
                break;
            case AggModels::PriPartList_ID:
                // Primary-particle list is known.  Use that
                // to construct the output sphere-tree.
                ppp = dynamic_cast<const AggModels::PriPartPrimary*>(sp.Primary());
                if (ppp != NULL) constructAgg_FM(*ppp);
                break;
        }
    }
}

// Constructs a PNode sphere-tree aggregate from the given 
// pri-part list primary using free-molecular collision dynamics.
void ParticleImage::constructAgg_FM(const AggModels::PriPartPrimary &pri)
{
    unsigned int i = 0;

    // Get number of primary particles.
    unsigned int n = pri.PriCount();

    // Create a vector of primary radii.
    fvector radii;
    for (unsigned int i=0; i!=n; ++i) {
        radii.push_back(pri.PriDiameter(i)*0.5e9); // Convert to nm.
    }

    // Randomly add the primaries to the image aggregate tree.
    m_root.Clear();
    while (radii.size() > 0) {
        unsigned int j = irnd(0, radii.size())-1;
        m_root.Insert(radii[j]);
        radii.erase(radii.begin()+j);
    }

    // Use the free-molecular regime to calculate the
    // aggregate structure.
    calc_FM(m_root);
}

// Constructs a PNode sphere-tree aggregate with uniform 
// primaries (equal diameter).  The diameter and primary
// count are passed as arguments.
void ParticleImage::uniformAgg_FM(unsigned int n, real d)
{
    real r = d * 0.5e9; // Convert to nm.

    // Add n identical primaries to the image data structure.
    m_root.Clear();
    for (unsigned int i=0; i!=n; ++i) {
        m_root.Insert(r);
    }

    // Use the free-molecular regime to calculate the
    // aggregate structure.
    calc_FM(m_root);
}

// Calculates the aggregate structure down from the given
// node.  Assumes that the tree leaves have been initialised
// with the correct radii, and recalculates their positions.
void ParticleImage::calc_FM(ImgNode &node)
{
    ImgNode *target = node.m_left;
    ImgNode *bullet = node.m_right;
    real sumrsqr = 0.0;

    if ((target != NULL) && (bullet != NULL)) {
        // Pass calculation down binary tree left & right branches.
        calc_FM(*target);
        calc_FM(*bullet);

        // The first part of the collision algorithm is to
        // randomly orientate both left and right aggregates.
        // They are both then placed so that their bounding
        // spheres are at the origin.

        // Rotate left node randomly about CoM.
        real phi1   = rnd() * 2.0 * PI;
        real theta1 = ((2.0*rnd())-1.0) * PI;
        target->RotateCOM(theta1, phi1);
        
        // Rotate right node randomly about CoM.
        real phi2   = rnd() * 2.0 * PI;
        real theta2 = ((2.0*rnd())-1.0) * PI;
        bullet->RotateCOM(theta2, phi2);

        // Move both spheres so that the bounding spheres
        // sit at the origin.
        target->CentreBoundSph();
        bullet->CentreBoundSph();

        // Perform the collision of the left and right nodes.
        // This ma require several iterations if the chosen
        // x-y displacement means that the aggregates cannot
        // collide in the z-direction.
        Coords::Vector D;
        real dxsqr=0.0, dysqr=0.0;
        real sumr=0.0, sumd=0.0;
        real dz1=0.0, dz2=0.0;
        bool hit = false;
        while (!hit) {
            // Need to reset target and bullet here, in case
            // the have been changed by the tree traversal
            // code below.
            target = node.m_left;
            bullet = node.m_right;
            sumr = target->Radius() + bullet->Radius();
            sumd = 2.0 * sumr;

            // Create a random displacement of the bullet node
            // in the x-y plane.  The squares are stored for
            // efficient computation.  The displacement is never
            // greater than the sum of the radii, therefore they
            // should always touch.
            D[0] = (rnd() * sumd) - 1.0; dxsqr = D[0]*D[0];
            D[1] = (rnd() * sumd) - 1.0; dysqr = D[1]*D[1];

            // Calculate the z-position for the collision of the
            // target and bullet.  We do this in case both the
            // target and bullet are leaf nodes, and hence the
            // binary tree traversal won't happen.
            hit = calcCollZ(sumr*sumr, dxsqr, dysqr, dz1);
            if (!hit) continue; // Should never happen.
            D[2] = target->m_cen_bsph[2] + dz1;

            // The next code determines the displacement along the z-axis
            // required for the target and bullet aggregates to touch.  This
            // requires falling down the tree progressively recalculating
            // the nearest nodes at each level, until the leaf nodes
            // are reached.

            while (hit) {
                // This next code calculates the minimum distance between the
                // target's children and bullet's children, or the target or
                // bullet if they have no children.  The two children with
                // the smallest separation are chosen as the next target
                // and bullet.
                if (target->IsLeaf()) {
                    if (bullet->IsLeaf()) {
                        // Check distance between leaf and bullet.
                        hit = calcCollZ(sumr*sumr, dxsqr, dysqr, dz1);
                        if (hit) D[2] = target->m_cen_bsph[2] + dz1;
                    } else {
                        // Check distance to bullet->left.
                        sumr = target->Radius() + bullet->m_left->Radius();
                        hit = calcCollZ(sumr*sumr, dxsqr, dysqr, dz1);
                        // Check distance to bullet->right.
                        sumr = target->Radius() + bullet->m_right->Radius();
                        hit = hit || calcCollZ(sumr*sumr, dxsqr, dysqr, dz2);
                        // Determine outcome.
                        if (hit) {
                            // Calculate minimum distance between aggregates.
                            D[2] = target->m_cen_bsph[2] + min(dz1, dz2);
                            // Choose next bullet.
                            if (dz1 <= dz2) {
                                bullet = bullet->m_left;
                            } else {
                                bullet = bullet->m_right;
                            }
                        }
                    }
                } else {
                    if (bullet->IsLeaf()) {
                        // Check distance to target->left.
                        sumr = target->m_left->Radius() + bullet->Radius();
                        hit = calcCollZ(sumr*sumr, dxsqr, dysqr, dz1);
                        // Check distance to bullet->right.
                        sumr = target->m_right->Radius() + bullet->Radius();
                        hit = hit||calcCollZ(sumr*sumr, dxsqr, dysqr, dz2);
                        // Determine outcome.
                        if (hit) {
                            // Calculate minimum distance between aggregates
                            // after choosing next target.
                            if (dz1 <= dz2) {
                                target = target->m_left;
                                D[2] = target->m_cen_bsph[2] + dz1;
                            } else {
                                target = target->m_right;
                                D[2] = target->m_cen_bsph[2] + dz2;
                            }
                        }
                    } else {
                        // Both the target and bullet have children.
                        // Need to check collision of all of them.  This code
                        // checks all combinations of child collisions
                        // consecutively, if the calculation results in a smaller
                        // separation then that is stored and the target/bullet
                        // combination is stored.
                        bool tleft=true, bleft=true;
                        // Calculate target->left and bullet->left.
                        sumr = target->m_left->Radius() + bullet->m_left->Radius();
                        hit = calcCollZ(sumr*sumr, dxsqr, dysqr, dz1);
                        // Calculate target->left and bullet->right.
                        sumr = target->m_left->Radius() + bullet->m_right->Radius();
                        hit = calcCollZ(sumr*sumr, dxsqr, dysqr, dz2);
                        if (dz2<dz1) {bleft=false; dz1=dz2;}
                        // Calculate target->right and bullet->left.
                        sumr = target->m_right->Radius() + bullet->m_left->Radius();
                        hit = calcCollZ(sumr*sumr, dxsqr, dysqr, dz2);
                        if (dz2<dz1) {tleft=false; bleft=true; dz1=dz2;}
                        // Calculate target->right and bullet->right.
                        sumr = target->m_right->Radius() + bullet->m_right->Radius();
                        hit = calcCollZ(sumr*sumr, dxsqr, dysqr, dz2);
                        if (dz2<dz1) {tleft=false; bleft=false; dz1=dz2;}
                        // Determine outcome.
                        if (hit) {
                            // Select next target and bullet.
                            if (tleft) {
                                target = target->m_left;
                            } else {
                                target = target->m_right;
                            }
                            if (bleft) {
                                bullet = bullet->m_left;
                            } else {
                                bullet = bullet->m_right;
                            }
                            // Calculate new bullet z-position.
                            D[2] = target->m_cen_bsph[2] + dz1;
                        }
                    }
                }
                if (target->IsLeaf() && bullet->IsLeaf()) break;
            } // while both target and bullet not leaves.
        } // While not hit.

        // We have a new location for the bullet (right node), so move it.
        node.m_right->Translate(D[0], D[1], D[2]);

        // Calculate properties of this node.
        node.CalcBoundSph();
        node.CalcCOM();
        node.CentreBoundSph();

    }
}

// Calculates the z-position of a bullet sphere for a +ve
// collision with a target sphere.  Returns true if the
// spheres collide, otherwise false.  A x-y displacement is
// included in the position of the bullet.
bool ParticleImage::calcCollZ(real sumrsqr, real dxsqr, real dysqr, real &dz)
{
    // Use Pythagoras theorem to calculate z-coord at which
    // spheres touch.
    dz = sumrsqr + dxsqr + dysqr; // Now contains dZ^2.
    if (dz >= 0.0) {
        // Spheres intersect.
        dz = sqrt(dz);
        return true;
    } else {
        // Spheres do not intersect.
        dz = 1.0e10; // A large number.
        return false;
    }
}
