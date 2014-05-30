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
	addPropertyChangeListener("Filter Type", this, &FilterDecimate_i::filterTypeChanged);
	addPropertyChangeListener("Number of Taps", this, &FilterDecimate_i::numberOfTapsChanged);
	addPropertyChangeListener("Center Frequency", this, &FilterDecimate_i::centerFrequencyChanged);
	addPropertyChangeListener("Bandwidth", this, &FilterDecimate_i::bandwidthChanged);
	addPropertyChangeListener("Gain", this, &FilterDecimate_i::gainChanged);
	addPropertyChangeListener("Output Rate", this, &FilterDecimate_i::outputRateChanged);

	//Initialize Filter Properties
	m_fc = 0.5;				// Filter cutoff frequency
	m_atten = 60.0f;   		// stop-band attenuation [dB]
	m_mu = 0.0f;       		// timing offset
	m_decimFactor = 1;		// decimation factor

	//Filter Objects Not Created
	m_complexFilterObjectCreated = false;
	m_realFilterObjectCreated = false;

	m_sriOut = bulkio::sri::create("BPSK_OUT"); //Create output sri object

	m_realFilter = NULL;
	m_complexFilter = NULL;
	_Input_Rate = 0;
	m_size = 0;
}

FilterDecimate_i::~FilterDecimate_i()
{
}

/************************* Property Change Listeners *******************************************/
void FilterDecimate_i::filterTypeChanged(const std::string *oldValue, const std::string *newValue)
{
	boost::mutex::scoped_lock lock(propertyLock_);

	//Only create filter if packet has been received
	if (_Input_Rate != 0){
		//Create filter if Filter Type changes
		if (_Data_Type == "complex") {
			createComplexFilter();
		} else if (_Data_Type == "real") {
			createRealFilter();
		}
	}
}

void FilterDecimate_i::numberOfTapsChanged(const unsigned short *oldValue, const unsigned short *newValue)
{
	boost::mutex::scoped_lock lock(propertyLock_);

	//Check that Number of Taps is valid
	if (Number_of_Taps < 1) {
		Number_of_Taps = *oldValue;
		std::cerr << "WARNING! -- 'Number_of_Taps' Must Be Greater Than 0!" << std::endl;
		std::cerr << "-- 'Number_of_Taps' set to " << *oldValue << std::endl;
	}else{
		Number_of_Taps = *newValue;
	}

	//Only create filter if packet has been received
	if (_Input_Rate != 0){
		//Create filter if Number of Taps change
		if (_Data_Type == "complex") {
			createComplexFilter();
		} else if (_Data_Type == "real") {
			createRealFilter();
		}
	}
}

void FilterDecimate_i::centerFrequencyChanged(const float *oldValue, const float *newValue)
{
	boost::mutex::scoped_lock lock(propertyLock_);

	//Check that Center_Frequency is valid
	if (Center_Frequency < 0) {
		Center_Frequency = *oldValue;
		std::cerr << "WARNING! -- 'Center_Frequency' Must Be Equal to or Greater Than 0!" << std::endl;
		std::cerr << "-- 'Center_Frequency' set to " << *oldValue << std::endl;
	}else{
		Center_Frequency = *newValue;
	}

	//Only create filter if packet has been received
	if (_Input_Rate !=0){
		//Create filter if Center Frequency changes
		if (_Data_Type == "complex") {
			createComplexFilter();
		} else if (_Data_Type == "real") {
			createRealFilter();
		}
	}
}

void FilterDecimate_i::bandwidthChanged(const float *oldValue, const float *newValue)
{
	boost::mutex::scoped_lock lock(propertyLock_);

	//Only create filter and set frequency cutoff if packet has been received
	if (_Input_Rate != 0){
		m_fc = Bandwidth / _Input_Rate;
		validateFilterCutoff();

		//Create filter if Bandwidth changes
		if (_Data_Type == "complex") {
			createComplexFilter();
		} else if (_Data_Type == "real") {
			createRealFilter();
		}
	}
}

void FilterDecimate_i::gainChanged(const float *oldValue, const float *newValue)
{
	boost::mutex::scoped_lock lock(propertyLock_);

	//Check that Gain is a valid
	if (Gain<=0) {
		Gain = *oldValue;
		std::cerr << "WARNING! -- 'Gain' Must Be Greater Than Zero!" << std::endl;
		std::cerr << "-- 'Gain' Set to " << *oldValue << std::endl;
	}else{
		Gain = *newValue;
	}

	//Only create filter if packet has been received
	if (_Input_Rate !=0){
		//Create filter if Gain changes
		if (_Data_Type == "complex") {
			createComplexFilter();
		} else if (_Data_Type == "real") {
			createRealFilter();
		}
	}
}

void FilterDecimate_i::outputRateChanged(const float *oldValue, const float *newValue)
{
	boost::mutex::scoped_lock lock(propertyLock_);

	//Only create filter and set decimation factor if packet has been received
	if (_Input_Rate != 0){

		//Set and validate new decimation factor
		m_decimFactor =  _Input_Rate / Output_Rate;
		validateDecimationFactor();

		//Resize output vector based on new decimation factor
		sizeVectors();

		//Create filter if outputRate changes
		if (_Data_Type == "complex") {
			createComplexFilter();
		} else if (_Data_Type == "real") {
			createRealFilter();
		}
	}
}

/***********************************************************************************************
************************************************************************************************/
int FilterDecimate_i::serviceFunction()
{
	//Input Data
	bulkio::InFloatPort::dataTransfer *input = dataFloat_in->getPacket(bulkio::Const::BLOCKING);

	if (not input) {return NOOP;}	// No data available

	{
		boost::mutex::scoped_lock lock(propertyLock_);

		//Determine if incoming data is real or complex
		if (input->SRI.mode) {
			_Data_Type = "complex";
			m_size = input->dataBuffer.size() / 2;
		} else {
			_Data_Type = "real";
			m_size = input->dataBuffer.size();
		}

		//Update variables and filters if input SRI changes
		if (input->sriChanged) {
			_Input_Rate = (1.0f / input->SRI.xdelta);

			//Calculate decimation factor and set filter cutoff
			m_decimFactor = _Input_Rate / Output_Rate;
			m_fc = Bandwidth / _Input_Rate;

			//Check that the values are valid
			validateDecimationFactor();
			validateFilterCutoff();

			//Create filter object
			if (_Data_Type == "complex") {
				createComplexFilter();
			} else if (_Data_Type == "real") {
				createRealFilter();
			}

			//Resize Vectors
			sizeVectors();

			//Update and push output SRI
			m_sriOut = input->SRI;
			m_sriOut.xdelta *= m_decimFactor;
			dataFloat_out->pushSRI(m_sriOut);
		}

		if (_Data_Type == "complex") {
			for (unsigned int i = 0; i < (m_size); ++i) {
				//Convert to c++ complex and filter
				m_inputComplex[i].real(input->dataBuffer[2 * i]);
				m_inputComplex[i].imag(input->dataBuffer[2 * i + 1]);
			}
			for (unsigned int i = 0; i < (m_size / m_decimFactor); i++) {
				//Filter complex data
				firdecim_crcf_execute(m_complexFilter, &m_inputComplex[m_decimFactor * i], &m_outputComplex[i], 0);

				//Convert Back to Interleaved complex, normalize output, multiply by gain
				m_output[2 * i] = real(m_outputComplex[i]) * Gain / m_decimFactor;
				m_output[2 * i + 1] = imag(m_outputComplex[i]) * Gain / m_decimFactor;
			}
		} else if (_Data_Type == "real") {

			for (unsigned int i = 0; i < (m_size / m_decimFactor); i++) {
				//Filter real data
				firdecim_rrrf_execute(m_realFilter, &input->dataBuffer[m_decimFactor * i], &m_output[i], 0);

				//normalize output & multiply by gain
				m_output[i] *= Gain / m_decimFactor;
			}
		}
	}

	//Push filtered and decimated data
	dataFloat_out->pushPacket(m_output, input->T, input->EOS, input->streamID);

	delete input;
	return NORMAL;
}

void FilterDecimate_i::createComplexFilter()
{
	//Destroy filter object if one already exists
	if (m_complexFilterObjectCreated) {firdecim_crcf_destroy(m_complexFilter);}
	float taps[Number_of_Taps];

	//Create filter coefficients
	liquid_firdes_kaiser(Number_of_Taps,m_fc,m_atten,m_mu,taps);

	//Shift filter taps by center frequency if bandpass filter
	if (Filter_Type == "bandpass") {
		std::complex<float> exp_arg;
		exp_arg.real(0);
		for(int i=0; i<Number_of_Taps; i++){
			exp_arg.imag(2.0*M_PI*Center_Frequency*i / _Input_Rate);
			taps[i] = real(taps[i]*exp(exp_arg));
		}
	}

	//Create filter object
	m_complexFilter = firdecim_crcf_create(m_decimFactor,taps,Number_of_Taps);
	m_complexFilterObjectCreated = true;
}

void FilterDecimate_i::createRealFilter()
{
	//Destroy filter object if one already exists
	if (m_realFilterObjectCreated) {firdecim_rrrf_destroy(m_realFilter);}
	float taps[Number_of_Taps];

	//Create filter coefficients
	liquid_firdes_kaiser(Number_of_Taps,m_fc,m_atten,m_mu,taps);

	if (Filter_Type == "bandpass") {
		float G = 1;
	    float fcLow = (Center_Frequency / _Input_Rate) - (m_fc/2.0);
	    float fcHigh = (Center_Frequency / _Input_Rate) + (m_fc/2.0);
		float T[Number_of_Taps];

	    int middleTap = (Number_of_Taps - 1) / 2;
	    double wcLow = 2 * M_PI * fcLow;
	    double wcHigh = 2 * M_PI * fcHigh;

	    //Shift filter taps by bandpass center frequency
	    for (int n = -middleTap; n <= middleTap; n++) {
			if (n == 0) {
				T[n + middleTap] = (wcHigh - wcLow) / M_PI * taps[n + middleTap];
			} else {
				T[n + middleTap] = (sin(n * wcHigh) - sin(n * wcLow)) / (n * M_PI) * taps[n + middleTap];
			}
		}

	    // Find the factor to normalize the gain, fmax.
	    // For band-pass, gain @ center freq = 1.0
	    double fmax = T[0 + middleTap];
	    for (int n = 1; n <= middleTap; n++) {
	    	fmax += 2 * T[n + middleTap] * cos(n * (wcLow + wcHigh) * 0.5);
	    }
	    G /= fmax;
	    // Normalize taps from calculated gain factor
	    for (int i = 0; i < Number_of_Taps; i++) {
	    	T[i] *= G;
	    }
	    memcpy(taps, T, Number_of_Taps*sizeof(float));
	}

	//Create filter object
	m_realFilter = firdecim_rrrf_create(m_decimFactor,taps,Number_of_Taps);
	m_realFilterObjectCreated = true;
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

void FilterDecimate_i::validateDecimationFactor()
{	//Decimation Factor is dependent on input SRI and Output_Rate property
	//Function allows filter cutoff to be checked without redundant code

	if (m_decimFactor < 1) {
		m_decimFactor = 1;
		Output_Rate = _Input_Rate;
		std::cerr << "WARNING! -- 'Output_Rate' Must Be Equal to or Smaller Than 'InputRate'!" << std::endl;
		std::cerr << "-- 'Output_Rate' set to 'InputRate'  --" << std::endl;
	}
}

void FilterDecimate_i::validateFilterCutoff()
{	//Filter Cutoff is dependent on input SRI and Bandwidth property
	//Function allows filter cutoff to be checked without redundant code

	if (m_fc > 0.5) {
		m_fc = 0.5;
		Bandwidth = 0.5 * _Input_Rate;
		std::cerr << "WARNING! -- 'Bandwidth' Must Be Equal to or Less Than ('InputRate' / 2)!" << std::endl;
		std::cerr << "-- 'Bandwidth' set to ('InputRate' / 2)  --" << std::endl;
	} else if (m_fc < 0) {
		m_fc = 0.5;
		Bandwidth = 0.5 * _Input_Rate;
		std::cerr << "WARNING! -- 'Bandwidth' Must Be Positive!" << std::endl;
		std::cerr << "-- 'Bandwidth' set to ('InputRate' / 2)  --" << std::endl;
	}
}

