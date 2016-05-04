#ifndef __FLUX_SRCS__
#define __FLUX_SRCS__

//! this fuction flags cells with active and passive flux sources (not to be confuesed with inter-cell numerical flux). It is initally called along with init_piles, then it is called everytime after adaptation is triggered by a begining of a flux-source.
extern void mark_flux_region(ElementsHashTable *ElmTable, NodeHashTable *NodeTable, MatProps *matprops, FluxProps *fluxprops,
                             TimeProps *timeprops);

//! this function triggers refinement when a flux sources starts adding material. TODO currently this is crudely implemented. Needs improvement, keith's idea of using binary flag is a good idea.
extern void adapt_fluxsrc_region(ElementsHashTable *ElemTable, NodeHashTable *NodeTable, MatProps *matprops,
                                 PileProps *pileprops, FluxProps *fluxprops, TimeProps *timeprops, double dt, int myid,
                                 int adaptflag);

// For Level Set Method (for future! when the level set method will be compatible with in flux sources)
//! this function triggers refinement when a flux sources starts adding material. TODO currently this is crudely implemented. Needs improvement, keith's idea of using binary flag is a good idea.
extern void adapt_fluxsrc_region_LevelSet(ElementsHashTable *ElemTable, NodeHashTable *NodeTable, MatProps *matprops,
                                 PileProps *pileprops, FluxProps *fluxprops, TimeProps *timeprops, double dt, int myid,
                                 int adaptflag);

//! this fuction calculates the flux contribution of the current cell based on its position relative to the source center and outflow-profie.
extern double calc_flux(Element *EmTemp, NodeHashTable *NodeTable, FluxProps *fluxprops, TimeProps *timeprops);

#endif
