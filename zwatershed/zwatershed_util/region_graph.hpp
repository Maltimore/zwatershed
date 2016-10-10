#pragma once

#include "types.hpp"

#include <cstddef>
#include <iostream>
#include <map>

template< typename ID, typename AG_P >
inline region_graph_ptr<ID,typename AG_P::element_type::element>
get_region_graph( const AG_P& aff_ptr,
                  const volume_ptr<ID> seg_ptr,
                  std::size_t max_segid)
{
    typedef typename AG_P::element_type AG;
    typedef typename AG::element        F;

    std::ptrdiff_t xdim = aff_ptr->shape()[1];
    std::ptrdiff_t ydim = aff_ptr->shape()[2];
    std::ptrdiff_t zdim = aff_ptr->shape()[3];

    volume<ID>& seg = *seg_ptr;
    auto& aff = *aff_ptr;

    region_graph_ptr<ID,F> rg_ptr( new region_graph<ID,F> );

    region_graph<ID,F>& rg = *rg_ptr;

    std::vector<std::map<ID,F>> edges(max_segid+1);

    for ( std::ptrdiff_t z = 0; z < zdim; ++z )
        for ( std::ptrdiff_t y = 0; y < ydim; ++y )
            for ( std::ptrdiff_t x = 0; x < xdim; ++x )
            {
                if ( (x > 0) && seg[x][y][z] && seg[x-1][y][z] )
                {
                    auto mm = std::minmax(seg[x][y][z], seg[x-1][y][z]);
                    F& curr = edges[mm.first][mm.second];
                    curr = std::max(curr, aff[0][x][y][z]);
                }
                if ( (y > 0) && seg[x][y][z] && seg[x][y-1][z] )
                {
                    auto mm = std::minmax(seg[x][y][z], seg[x][y-1][z]);
                    F& curr = edges[mm.first][mm.second];
                    curr = std::max(curr, aff[1][x][y][z]);
                }
                if ( (z > 0) && seg[x][y][z] && seg[x][y][z-1] )
                {
                    auto mm = std::minmax(seg[x][y][z], seg[x][y][z-1]);
                    F& curr = edges[mm.first][mm.second];
                    curr = std::max(curr, aff[2][x][y][z]);
                }
            }

    for ( ID id1 = 1; id1 <= max_segid; ++id1 )
    {
        for ( const auto& p: edges[id1] )
        {
            rg.emplace_back(p.second, id1, p.first);
            //std::cout << p.second << " " << id1 << " " << p.first << "\n";
        }
    }

    std::cout << "Region graph size: " << rg.size() << std::endl;

    std::stable_sort(std::begin(rg), std::end(rg),
                     std::greater<std::tuple<F,ID,ID>>());

    return rg_ptr;
}
