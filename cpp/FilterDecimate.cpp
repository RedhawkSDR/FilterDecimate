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

/**************************************************************************

    This is the component code. This file contains the child class where
    custom functionality can be added to the component. Custom
    functionality to the base class can be extended here. Access to
    the ports can also be done from this class

 	Source: FilterDecimate.spd.xml
 	Generated on: Thu Sep 26 16:01:49 EDT 2013
 	REDHAWK IDE
 	Version: 1.8.4
 	Build id: R201305151907

**************************************************************************/

#include "FilterDecimate.h"

PREPARE_LOGGING(FilterDecimate_i)

FilterDecimate_i::FilterDecimate_i(const char *uuid, const char *label) :
    FilterDecimate_base(uuid, label)
{
	//Enable property change listeners to call propertyChangeListener when a property changes
	setPropertyChangeListener("Filter Type", this,
			&FilterDecimate_i::propertyChangeListener);
	setPropertyChangeListener("Number of Taps", this,
			&FilterDecimate_i::propertyChangeListener);
	setPropertyChangeListener("Center Frequency", this,
			&FilterDecimate_i::propertyChangeListener);
	setPropertyChangeListener("Bandwidth", this,
			&FilterDecimate_i::propertyChangeListener);
	setPropertyChangeListener("Gain", this,
			&FilterDecimate_i::propertyChangeListener);
	setPropertyChangeListener("Output Rate", this,
			&FilterDecimate_i::propertyChangeListener);

	//Initialize Filter Properties
	m_fc = 0.5;				// filter cutoff
	m_atten = 60.0f;   		// stop-band attenuation [dB]
	m_mu = 0.0f;       		// timing offset
	m_decimFactor = 1;		// decimation factor

	//Filter Objects Not Created
	m_complexFilterObjectCreated = 0;
	m_realFilterObjectCreated = 0;
	m_propertyChanged = true;

	m_realFilter = 0;
	m_complexFilter = 0;
	m_delta = 0;
	m_size = 0;
}

FilterDecimate_i::~FilterDecimate_i()
{
}


/***********************************************************************************************
************************************************************************************************/
int FilterDecimate_i::serviceFunction()
{
	//Input Data
	bulkio::InFloatPort::dataTransfer *input = dataFloat_in->getPacket(bulkio::Const::BLOCKING);

	if (not input) {return NOOP;}	// No data available

	//Determine if incoming data is real or complex
	if (input->SRI.mode) {
		_Data_Type = "complex";
		m_size = input->dataBuffer.size()/2;
	} else {
		_Data_Type = "real";
		m_size = input->dataBuffer.size();
	}

	//Update variables and filters if sri or properties change
	if (input->sriChanged || m_propertyChanged) {
		m_delta = input->SRI.xdelta;
		_Input_Rate = (1.0f / m_delta);

		//Check that new properties are valid
		checkProperties();

		//Create filter if properties or SRI changed
		if (_Data_Type == "complex") {
			createComplexFilter();
		} else if (_Data_Type == "real") {
			createRealFilter();
		}

		sizeVectors();

		//Update sample rate and push SRI
		input->SRI.xdelta *= m_decimFactor;
		dataFloat_out->pushSRI(input->SRI);

		m_propertyChanged = false;
	}

	if (_Data_Type == "complex") {
		//Filter & Decimate Data
		for (unsigned int i=0; i< (m_size); ++i){
			//Convert to c++ complex and filter
			m_inputComplex[i].real(input->dataBuffer[2*i]);
			m_inputComplex[i].imag(input->dataBuffer[2*i+1]);
		}
		for (unsigned int i=0; i< (m_size/m_decimFactor); i++) {
		   	firdecim_crcf_execute(m_complexFilter, &m_inputComplex[m_decimFactor * i], &m_outputComplex[i], 0);
		   	//Convert Back to Interleaved complex and normalize output
		   	m_output[2*i] = real(m_outputComplex[i]) * Gain / m_decimFactor;
		   	m_output[2*i+1] = imag(m_outputComplex[i]) * Gain /m_decimFactor;
		}
	} else if (_Data_Type == "real") {
		//Filter & Decimate Data
		for (unsigned int i = 0; i < (m_size/m_decimFactor); i++) {
			firdecim_rrrf_execute(m_realFilter, &input->dataBuffer[m_decimFactor * i], &m_output[i], 0);
			m_output[i] *= Gain / m_decimFactor;	//normalize output & multiply by gain
		}
	}

	dataFloat_out->pushPacket(m_output, input->T, input->EOS, input->streamID);

	delete input;
	return NORMAL;
}

void FilterDecimate_i::createComplexFilter()
{
	//Destroy filter object if one already exists
	if (m_complexFilterObjectCreated) {firdecim_crcf_destroy(m_complexFilter);}

	float h[Number_of_Taps];		//Vector of taps
	m_decimFactor = 1/(Output_Rate*m_delta);	//Calculate decimation factor
	m_fc = m_delta * Bandwidth;		//Set filter cutoff

	//Make sure m_fc is valid
	if (m_fc > 0.5) {
		m_fc = 0.5;
		std::cout << "m_fc > 0.5, setting m_fc to 0.5" << std::endl;
	} else if (m_fc < 0) {
		m_fc = 0;
		std::cout << "m_fc < 0, setting m_fc to 0.0" << std::endl;
	}

	//Make sure decimation factor is valid
	if (m_decimFactor < 1) {
		m_decimFactor = 1;
		Output_Rate = 1/m_delta;
		std::cerr << "WARNING! -- 'Output_Rate' Must Be Equal to or Smaller Than 'InputRate'!" << std::endl;
		std::cerr << "-- 'Output_Rate' set to 'InputRate'  --" << std::endl;
	}

	//Create filter coefficients
	liquid_firdes_kaiser(Number_of_Taps,m_fc,m_atten,m_mu,h);

	//Shift filter taps by center frequency if bandpass filter
	if (Filter_Type == "bandpass") {
		std::complex<float> exp_arg;
		exp_arg.real(0);
		for(int i=0; i<Number_of_Taps; i++){
			exp_arg.imag(2.0*M_PI*Center_Frequency*i*m_delta);
			h[i] = real(h[i]*exp(exp_arg));
		}
	}

	//Create filter object
	m_complexFilter = firdecim_crcf_create(m_decimFactor,h,Number_of_Taps);
	m_complexFilterObjectCreated = 1;
}

void FilterDecimate_i::createRealFilter()
{
	//Destroy filter object if one already exists
	if (m_realFilterObjectCreated) {firdecim_rrrf_destroy(m_realFilter);}
	float h[Number_of_Taps];		//Vector of taps
	m_decimFactor = 1/(Output_Rate*m_delta);	//Calculate decimation factor
	m_fc = m_delta * Bandwidth;		//Set filter cutoff

	Output_Rate = 1/(m_decimFactor*m_delta);

	//Make sure m_fc is valid
	if (m_fc > 0.5) {
		m_fc = 0.5;
		std::cout << "m_fc > 0.5, setting m_fc to 0.5" << std::endl;
	} else if (m_fc < 0) {
		m_fc = 0;
		std::cout << "m_fc < 0, setting m_fc to 0.0" << std::endl;
	}

	//Make sure decimation factor is valid
	if (m_decimFactor < 1) {
		m_decimFactor = 1;
		Output_Rate = 1/m_delta;
		std::cerr << "WARNING! -- 'Output_Rate' Must Be Equal to or Smaller Than 'InputRate'!" << std::endl;
		std::cerr << "-- 'Output_Rate' set to 'InputRate'  --" << std::endl;
	}

	//Create filter coefficients
	liquid_firdes_kaiser(Number_of_Taps,m_fc,m_atten,m_mu,h);

	if (Filter_Type == "bandpass") {
		float G = 1;
	    float fcLow = Center_Frequency*m_delta - (m_fc/2.0);
	    float fcHigh = Center_Frequency*m_delta + (m_fc/2.0);
		float T[Number_of_Taps];

	    int M = (Number_of_Taps - 1) / 2;
	    double fwT0 = 2 * M_PI * fcLow;
	    double fwT1 = 2 * M_PI * fcHigh;

	    //Shift filter taps by bandpass center frequency
	    for (int n = -M; n <= M; n++) {
			if (n == 0) {
				T[n + M] = (fwT1 - fwT0) / M_PI * h[n + M];
			} else {
				T[n + M] = (sin(n * fwT1) - sin(n * fwT0)) / (n * M_PI) * h[n+ M];
			}
		}

	    // Find the factor to normalize the gain, fmax.
	    // For band-pass, gain @ center freq = 1.0
	    double fmax = T[0 + M];
	    for (int n = 1; n <= M; n++) {
	    	fmax += 2 * T[n + M] * cos(n * (fwT0 + fwT1) * 0.5);
	    }
	    G /= fmax;
	    // Normalize taps from calculated gain factor
	    for (int i = 0; i < Number_of_Taps; i++) {
	    	T[i] *= G;
	    }
	    memcpy(h, T, Number_of_Taps*sizeof(float));
	}

	//Create filter object
	m_realFilter = firdecim_rrrf_create(m_decimFactor,h,Number_of_Taps);
	m_realFilterObjectCreated = 1;
}

void FilterDecimate_i::sizeVectors()
{
	if (_Data_Type == "complex") {
		//Resize vectors for complex data
		m_inputComplex.resize(m_size);
		m_outputComplex.resize(m_size/m_decimFactor);
		m_output.resize(m_outputComplex.size() * 2);
	}
	if (_Data_Type == "real") {
		//Resize vectors for complex data
		m_output.resize(m_size/m_decimFactor);
	}
}

void FilterDecimate_i::propertyChangeListener(const std::string& id)
{
	m_propertyChanged = true;
}

void FilterDecimate_i::checkProperties()
{
	//Check that Gain is a valid
	if (Gain<=0) {
		Gain = 1;
		std::cerr << "WARNING! -- 'Gain' Must Be Greater Than Zero!" << std::endl;
		std::cerr << "-- 'Gain' Set to 1" << std::endl;
	}

	//Check that Bandwidth is valid
	if (Bandwidth > (Output_Rate/2)) {
		Bandwidth = Output_Rate/2;
		std::cerr << "WARNING! -- 'Bandwidth' Must Be Equal to or Less Than 0.5*'OuputRate'!" << std::endl;
		std::cerr << "-- 'Bandwidth' Set to 0.5*'Output_Rate'" << std::endl;
	}

	//Check that Center_Frequency is valid
	if (Center_Frequency < 0) {
		Center_Frequency = 1;
		Filter_Type == "lowpass";
		std::cerr << "WARNING! -- 'Center_Frequency' Must Be Equal to or Greater Than 0!" << std::endl;
		std::cerr << "-- 'Center_Frequency' set to 0 && 'Filter_Type' set to 'lowpass'" << std::endl;
	}

	//Check that Taps is valid
	if (Number_of_Taps < 1) {
		Number_of_Taps = 20;
		std::cerr << "WARNING! -- 'Number_of_Taps' Must Be Greater Than 0!" << std::endl;
		std::cerr << "-- 'Number_of_Taps' set to 20" << std::endl;
	}
}

