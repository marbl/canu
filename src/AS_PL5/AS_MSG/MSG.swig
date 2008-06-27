/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author: Guillaume Marcais <guillaume@marcais.net>
 *************************************************************************/

%module "AS::MSG"
%{
#include "AS_MSG/AS_MSG_pmesg.h"
#include "MSG.h"
%}

/* AS pmsg definitions
 */
%include pmesg.h

/* Parser definitions
 */
%include parser.h

/* 3 letters acronymes of struct members
 * TODO: This should really be moved to Parser.pm.
 */
%include acronyms.h
