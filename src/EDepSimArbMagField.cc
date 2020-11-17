// License and Disclaimer
//
// For use with software authored by the Geant4 Collaboration
//
// MIT License
//
// Copyright (c) 2020 Andrew Cudd
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
//
// Class for storing and handling an arbitrary magnetic field.
//
// History:
// - 2020.04.14 A.Cudd created
// - 2020.07.28 C.McGrew updated license with permission of A.Cudd
// - 2020.11.17 A.Cudd added code to read in the units for the field
//              position and strength (e.g. cm and tesla) along with a
//              flag to indicate if the map is symmetric about the XYZ
//              axes. More in the README.
//
// -------------------------------------------------------------------

#include "EDepSimArbMagField.hh"

EDepSim::ArbMagField::ArbMagField()
{
}

bool EDepSim::ArbMagField::ReadFile(const std::string& fname)
{
    m_filename = fname;
    m_is_symmetric = false;
    std::fstream fin(fname, std::fstream::in);

    if(!fin.is_open())
    {
        EDepSimError("Can't read " << fname << std::endl);
        return false;
    }
    else
    {
        EDepSimLog("Reading " << fname << " ...");
        int xcount{-1}, ycount{-1}, zcount{-1};
        double xcurr{0}, ycurr{0}, zcurr{0};
        std::string line;

        std::string pos_units("mm");
        std::string field_units("tesla");

        while(std::getline(fin >> std::ws, line))
        {
            std::stringstream ss(line);

            if(ss.str().front() == '#')
                continue;

            ss >> m_offset[0] >> m_offset[1] >> m_offset[2]
                >> m_delta[0] >> m_delta[1] >> m_delta[2]
                >> pos_units >> field_units
                >> std::boolalpha >> m_is_symmetric;
            break;
        }

        m_position_units = G4UnitDefinition::GetValueOf(pos_units);
        m_field_units = G4UnitDefinition::GetValueOf(field_units);

        for(auto& val : m_offset)
            val *= m_position_units;

        for(auto& val : m_delta)
            val *= m_position_units;

        while(std::getline(fin >> std::ws, line))
        {
            double x{0}, y{0}, z{0}, fx{0}, fy{0}, fz{0}, f{0};
            std::stringstream ss(line);

            if(ss.str().front() == '#')
                continue;

            ss >> x >> y >> z >> fx >> fy >> fz >> f;

            x *= m_position_units;
            y *= m_position_units;
            z *= m_position_units;

            if(std::abs(x - xcurr) > 0.0 || xcount < 0)
            {
                xcurr = x;
                xcount += 1;
                ycount = -1;
                m_field_x.emplace_back(std::vector<std::vector<double>>{});
                m_field_y.emplace_back(std::vector<std::vector<double>>{});
                m_field_z.emplace_back(std::vector<std::vector<double>>{});
            }

            if(std::abs(y - ycurr) > 0.0 || ycount < 0)
            {
                ycurr = y;
                ycount += 1;
                zcount = -1;
                m_field_x[xcount].emplace_back(std::vector<double>{});
                m_field_y[xcount].emplace_back(std::vector<double>{});
                m_field_z[xcount].emplace_back(std::vector<double>{});
            }

            m_field_x[xcount][ycount].push_back(fx * m_field_units);
            m_field_y[xcount][ycount].push_back(fy * m_field_units);
            m_field_z[xcount][ycount].push_back(fz * m_field_units);

            if(std::abs(z - zcurr) > 0.0 || zcount < 0)
            {
                zcurr = z;
                zcount += 1;
            }
        }
    }

    return true;
}

void EDepSim::ArbMagField::GetFieldValue(const G4double pos[4], G4double* field) const
{
    const double x = m_is_symmetric ? std::abs(pos[0]) : pos[0];
    const double y = m_is_symmetric ? std::abs(pos[1]) : pos[1];
    const double z = m_is_symmetric ? std::abs(pos[2]) : pos[2];

    EDepSim::Cubic interp;
    field[0] = interp.interpolate(x, y, z, m_field_x, m_delta[0], m_delta[1], m_delta[2], m_offset[0], m_offset[1], m_offset[2]);
    field[1] = interp.interpolate(x, y, z, m_field_y, m_delta[0], m_delta[1], m_delta[2], m_offset[0], m_offset[1], m_offset[2]);
    field[2] = interp.interpolate(x, y, z, m_field_z, m_delta[0], m_delta[1], m_delta[2], m_offset[0], m_offset[1], m_offset[2]);
}

void EDepSim::ArbMagField::PrintInfo() const
{
    EDepSimLog("Printing values for magnetic field.");
    EDepSimLog("m_filename : " << m_filename
              << "\nm_offset   : " << m_offset[0] << ", " << m_offset[1] << ", " << m_offset[2]
              << "\nm_delta    : " << m_delta[0] << ", " << m_delta[1] << ", " << m_delta[2]
              << "\nm_symmetric: " << std::boolalpha << m_is_symmetric);
}
