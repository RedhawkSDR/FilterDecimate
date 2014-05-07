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

#include "FilterDecimate_base.h"

/*******************************************************************************************

    AUTO-GENERATED CODE. DO NOT MODIFY

    The following class functions are for the base class for the component class. To
    customize any of these functions, do not modify them here. Instead, overload them
    on the child class

******************************************************************************************/

FilterDecimate_base::FilterDecimate_base(const char *uuid, const char *label) :
    Resource_impl(uuid, label),
    ThreadedComponent()
{
    loadProperties();

    dataFloat_in = new bulkio::InFloatPort("dataFloat_in");
    addPort("dataFloat_in", dataFloat_in);
    dataFloat_out = new bulkio::OutFloatPort("dataFloat_out");
    addPort("dataFloat_out", dataFloat_out);
}

FilterDecimate_base::~FilterDecimate_base()
{
    delete dataFloat_in;
    dataFloat_in = 0;
    delete dataFloat_out;
    dataFloat_out = 0;
}

/*******************************************************************************************
    Framework-level functions
    These functions are generally called by the framework to perform housekeeping.
*******************************************************************************************/
void FilterDecimate_base::start() throw (CORBA::SystemException, CF::Resource::StartError)
{
    Resource_impl::start();
    ThreadedComponent::startThread();
}

void FilterDecimate_base::stop() throw (CORBA::SystemException, CF::Resource::StopError)
{
    Resource_impl::stop();
    if (!ThreadedComponent::stopThread()) {
        throw CF::Resource::StopError(CF::CF_NOTSET, "Processing thread did not die");
    }
}

void FilterDecimate_base::releaseObject() throw (CORBA::SystemException, CF::LifeCycle::ReleaseError)
{
    // This function clears the component running condition so main shuts down everything
    try {
        stop();
    } catch (CF::Resource::StopError& ex) {
        // TODO - this should probably be logged instead of ignored
    }

    Resource_impl::releaseObject();
}

void FilterDecimate_base::loadProperties()
{
    addProperty(_Data_Type,
                "_Data Type",
                "_Data Type",
                "readonly",
                "",
                "external",
                "configure");

    addProperty(_Input_Rate,
                "_Input Rate",
                "_Input Rate",
                "readonly",
                "S/s",
                "external",
                "configure");

    addProperty(Bandwidth,
                100000,
                "Bandwidth",
                "",
                "readwrite",
                "Hz",
                "external",
                "configure");

    addProperty(Center_Frequency,
                0,
                "Center Frequency",
                "",
                "readwrite",
                "Hz",
                "external",
                "configure");

    addProperty(Filter_Type,
                "lowpass",
                "Filter Type",
                "",
                "readwrite",
                "",
                "external",
                "configure");

    addProperty(Gain,
                1,
                "Gain",
                "",
                "readwrite",
                "",
                "external",
                "configure");

    addProperty(Number_of_Taps,
                20,
                "Number of Taps",
                "",
                "readwrite",
                "",
                "external",
                "configure");

    addProperty(Output_Rate,
                256000,
                "Output Rate",
                "",
                "readwrite",
                "S/s",
                "external",
                "configure");

}


