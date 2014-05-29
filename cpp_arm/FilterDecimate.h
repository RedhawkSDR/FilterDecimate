/*
 * This file is protected by Copyright. Please refer to the COPYRIGHT file
 * distributed with this source distribution.
 *
 * This file is part of FilterDecimate.
 *
 * FilterDecimate is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * FilterDecimate is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/.
 */

#ifndef FILTERDECIMATE_IMPL_H
#define FILTERDECIMATE_IMPL_H

#include "FilterDecimate_base.h"
#include <liquid.h>
#include <complex>

class FilterDecimate_i;

class FilterDecimate_i : public FilterDecimate_base
{
    ENABLE_LOGGING
    public:
        FilterDecimate_i(const char *uuid, const char *label);
        ~FilterDecimate_i();
        int serviceFunction();
    private:
        //Data Vectors
        std::vector< std::complex<float> > m_inputComplex;
        std::vector< std::complex<float> > m_outputComplex;
		std::vector<float> m_output;

		//Filter Objects
		firdecim_crcf m_complexFilter;
		firdecim_rrrf m_realFilter;
		int m_complexFilterObjectCreated;
		int m_realFilterObjectCreated;
		bool m_propertyChanged;

		//Filter Properties
		float m_fc;				// filter cutoff
		float m_atten;         	// stop-band attenuation [dB]
		float m_mu;          	// timing offset
		int m_decimFactor;		// decimation factor

		//SRI Data
		double m_delta;
		float m_size;

		//Member Functions
		void createComplexFilter();
		void createRealFilter();
		void checkProperties();
		void sizeVectors();
		void propertyChangeListener(const std::string&);
};

#endif
