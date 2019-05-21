module ExternalModelPFLOTRANMod
#ifdef USE_PETSC_LIB

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscviewer.h"

  use abortutils                   , only : endrun
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use shr_log_mod                  , only : errMsg => shr_log_errMsg
  use EMI_DataMod                  , only : emi_data_list, emi_data
  use clm_varctl                   , only : iulog
  use ExternalModelBaseType        , only : em_base_type
  use MultiPhysicsProbVSFM         , only : mpp_vsfm_type
  use decompMod                    , only : bounds_type
  use pflotran_model_module
  use ExternalModelConstants
  use EMI_Atm2LndType_Constants
  use EMI_CanopyStateType_Constants
  use EMI_ColumnType_Constants
  use EMI_EnergyFluxType_Constants
  use EMI_Filter_Constants
  use EMI_Landunit_Constants
  use EMI_SoilHydrologyType_Constants
  use EMI_SoilStateType_Constants
  use EMI_TemperatureType_Constants
  use EMI_WaterFluxType_Constants
  use EMI_WaterStateType_Constants
  use clm_pflotran_interface_data
  !
  implicit none
  !

  type, public, extends(em_base_type) :: em_pflotran_type
     ! ----------------------------------------------------------------------
     ! Indicies required during the initialization
     ! ----------------------------------------------------------------------
     integer :: index_l2e_init_col_active
     integer :: index_l2e_init_col_type
     integer :: index_l2e_init_col_landunit_index
     integer :: index_l2e_init_col_gridcell_index
     integer :: index_l2e_init_col_zi
     integer :: index_l2e_init_col_dz
     integer :: index_l2e_init_col_z
     integer :: index_l2e_init_col_area

     integer :: index_l2e_init_state_wtd
     integer :: index_l2e_init_state_soilp

     integer :: index_l2e_init_h2osoi_liq
     integer :: index_l2e_init_h2osoi_ice

     integer :: index_e2l_init_state_h2osoi_liq
     integer :: index_e2l_init_state_h2osoi_ice
     integer :: index_e2l_init_state_h2osoi_vol
     integer :: index_e2l_init_state_wtd
     integer :: index_e2l_init_parameter_watsatc
     integer :: index_e2l_init_parameter_hksatc
     integer :: index_e2l_init_parameter_bswc
     integer :: index_e2l_init_parameter_sucsatc

     integer :: index_e2l_init_flux_mflx_snowlyr_col
     integer :: index_l2e_init_flux_mflx_snowlyr_col

     integer :: index_l2e_init_landunit_type
     integer :: index_l2e_init_landunit_lakepoint
     integer :: index_l2e_init_landunit_urbanpoint

     integer :: index_l2e_init_parameter_watsatc
     integer :: index_l2e_init_parameter_hksatc
     integer :: index_l2e_init_parameter_bswc
     integer :: index_l2e_init_parameter_sucsatc
     integer :: index_l2e_init_parameter_effporosityc

     ! ----------------------------------------------------------------------
     ! Indicies required during timestepping
     ! ----------------------------------------------------------------------

     integer :: index_l2e_state_tsoil
     integer :: index_l2e_state_h2osoi_liq
     integer :: index_l2e_state_h2osoi_ice

     integer :: index_e2l_state_h2osoi_liq
     integer :: index_e2l_state_h2osoi_ice
     integer :: index_e2l_state_wtd
     integer :: index_e2l_state_soilp

     integer :: index_l2e_flux_infil
     integer :: index_l2e_flux_et
     integer :: index_l2e_flux_dew
     integer :: index_l2e_flux_snow_sub
     integer :: index_l2e_flux_snowlyr
     integer :: index_l2e_flux_drainage

     integer :: index_e2l_flux_qrecharge

     integer :: index_l2e_filter_hydrologyc
     integer :: index_l2e_filter_num_hydrologyc

     integer :: index_l2e_column_zi

     ! IDs to indentify the conditions for VSFM
     integer :: vsfm_cond_id_for_infil
     integer :: vsfm_cond_id_for_et
     integer :: vsfm_cond_id_for_dew
     integer :: vsfm_cond_id_for_drainage
     integer :: vsfm_cond_id_for_snow
     integer :: vsfm_cond_id_for_sublimation
     integer :: vsfm_cond_id_for_lateral_flux

     type(pflotran_model_type), pointer :: pflotran_m

   contains

     procedure, public :: Populate_L2E_Init_List  => EM_PFLOTRAN_Populate_L2E_Init_List
     procedure, public :: Populate_E2L_Init_List  => EM_PFLOTRAN_Populate_E2L_Init_List
     procedure, public :: Populate_L2E_List       => EM_PFLOTRAN_Populate_L2E_List
     procedure, public :: Populate_E2L_List       => EM_PFLOTRAN_Populate_E2L_List
     procedure, public :: PreInit                 => EM_PFLOTRAN_PreInit
     procedure, public :: Init                    => EM_PFLOTRAN_Init
     procedure, public :: Solve                   => EM_PFLOTRAN_Solve

  end type em_pflotran_type

contains

  !------------------------------------------------------------------------
  subroutine EM_PFLOTRAN_Populate_L2E_Init_List(this, l2e_init_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by PFLOTRAN from ALM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_pflotran_type)             :: this
    class(emi_data_list), intent(inout) :: l2e_init_list
    !
    ! !LOCAL VARIABLES:
    integer              , pointer       :: em_stages(:)
    integer                              :: number_em_stages
    integer                              :: id
    integer                              :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_INITIALIZATION_STAGE

    id                                         = L2E_STATE_WTD
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_wtd              = index

    id                                         = L2E_STATE_VSFM_PROGNOSTIC_SOILP
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_soilp            = index

    id                                         = L2E_FLUX_RESTART_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_flux_mflx_snowlyr_col  = index

    id                                         = L2E_COLUMN_ACTIVE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_active             = index

    id                                         = L2E_COLUMN_TYPE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_type               = index

    id                                         = L2E_COLUMN_LANDUNIT_INDEX
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_landunit_index     = index

    id                                         = L2E_COLUMN_GRIDCELL_INDEX
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_gridcell_index     = index

    id                                         = L2E_COLUMN_ZI
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_zi                 = index

    id                                         = L2E_COLUMN_DZ
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_dz                 = index

    id                                         = L2E_COLUMN_Z
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_z                  = index

    id                                         = L2E_COLUMN_AREA
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_area               = index

    id                                         = L2E_LANDUNIT_TYPE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_type          = index

    id                                         = L2E_LANDUNIT_LAKEPOINT
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_lakepoint     = index

    id                                         = L2E_LANDUNIT_URBANPOINT
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_urbanpoint    = index

    id                                         = L2E_PARAMETER_WATSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_watsatc      = index

    id                                         = L2E_PARAMETER_HKSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_hksatc       = index

    id                                         = L2E_PARAMETER_BSWC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_bswc         = index

    id                                         = L2E_PARAMETER_SUCSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_sucsatc      = index

    id                                         = L2E_PARAMETER_EFFPOROSITYC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_effporosityc = index

    id                                        = L2E_STATE_H2OSOI_LIQ_NLEVGRND
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_h2osoi_liq            = index

    id                                        = L2E_STATE_H2OSOI_ICE_NLEVGRND
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_h2osoi_ice            = index

    deallocate(em_stages)

  end subroutine EM_PFLOTRAN_Populate_L2E_Init_List

  !------------------------------------------------------------------------
  subroutine EM_PFLOTRAN_Populate_E2L_Init_List(this, e2l_init_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by PFLOTRAN from ALM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_pflotran_type)             :: this
    class(emi_data_list), intent(inout) :: e2l_init_list
    !
    ! !LOCAL VARIABLES:
    integer              , pointer       :: em_stages(:)
    integer                              :: number_em_stages
    integer                              :: id
    integer                              :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_INITIALIZATION_STAGE

    id                                        = E2L_STATE_H2OSOI_LIQ
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_liq      = index

    id                                        = E2L_STATE_H2OSOI_ICE
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_ice      = index

    id                                        = E2L_STATE_H2OSOI_VOL_NLEVGRND
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_vol      = index

    id                                        = E2L_STATE_WTD
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_wtd             = index

    id                                        = E2L_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_flux_mflx_snowlyr_col = index

    id                                         = E2L_PARAMETER_WATSATC
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_parameter_watsatc      = index

    id                                         = E2L_PARAMETER_HKSATC
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_parameter_hksatc       = index

    id                                         = E2L_PARAMETER_BSWC
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_parameter_bswc         = index

    id                                         = E2L_PARAMETER_SUCSATC
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_parameter_sucsatc      = index

    deallocate(em_stages)

  end subroutine EM_PFLOTRAN_Populate_E2L_Init_List

  !------------------------------------------------------------------------
  subroutine EM_PFLOTRAN_Populate_L2E_List(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by PFLOTRAN from ALM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_pflotran_type)             :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !
    !
    ! !LOCAL VARIABLES:
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_PFLOTRAN_SOIL_HYDRO_STAGE

    id                                   = L2E_STATE_TSOIL_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_tsoil           = index

    id                                   = L2E_STATE_H2OSOI_LIQ_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_liq      = index

    id                                   = L2E_STATE_H2OSOI_ICE_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_ice      = index

    id                                   = L2E_FLUX_INFIL_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_infil            = index

    id                                   = L2E_FLUX_VERTICAL_ET_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_et               = index

    id                                   = L2E_FLUX_DEW_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_dew              = index

    id                                   = L2E_FLUX_SNOW_SUBLIMATION_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_snow_sub         = index

    id                                   = L2E_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_snowlyr          = index

    id                                   = L2E_FLUX_DRAINAGE_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_drainage         = index

    id                                   = L2E_FILTER_HYDROLOGYC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_hydrologyc     = index

    id                                   = L2E_FILTER_NUM_HYDROLOGYC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_num_hydrologyc = index

    id                                   = L2E_COLUMN_ZI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_column_zi             = index

    deallocate(em_stages)

  end subroutine EM_PFLOTRAN_Populate_L2E_List

  !------------------------------------------------------------------------
  subroutine EM_PFLOTRAN_Populate_E2L_List(this, e2l_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by PFLOTRAN from ALM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_pflotran_type)             :: this
    class(emi_data_list), intent(inout) :: e2l_list
    !
    ! !LOCAL VARIABLES:
    integer              , pointer       :: em_stages(:)
    integer                              :: number_em_stages
    integer                              :: id
    integer                              :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_PFLOTRAN_SOIL_HYDRO_STAGE

    id                              = E2L_STATE_H2OSOI_LIQ
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_h2osoi_liq = index

    id                              = E2L_STATE_H2OSOI_ICE
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_h2osoi_ice = index

    id                              = E2L_STATE_WTD
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_wtd        = index

    id                              = E2L_STATE_VSFM_PROGNOSTIC_SOILP
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_soilp      = index

    id                              = E2L_FLUX_AQUIFER_RECHARGE
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_qrecharge   = index

    deallocate(em_stages)

  end subroutine EM_PFLOTRAN_Populate_E2L_List

  !------------------------------------------------------------------------
  subroutine EM_PFLOTRAN_PreInit(this, prefix, restart_stamp)
    !
    use spmdMod     , only : mpicom
    !
    implicit none
    !
    class(em_pflotran_type)      :: this
    character(len=*), intent(in) :: prefix
    character(len=*), intent(in) :: restart_stamp
    
    this%pflotran_m => pflotranModelCreate(mpicom, prefix)
    call pflotranModelSetupRestart(this%pflotran_m, restart_stamp)

  end subroutine EM_PFLOTRAN_PreInit

  !------------------------------------------------------------------------
  subroutine EM_PFLOTRAN_Init(this, l2e_init_list, e2l_init_list, iam, bounds_clump)

    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_pflotran_type)              :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    integer              , intent(in)    :: iam
    type(bounds_type)    , intent(in)    :: bounds_clump
    !
    ! LOCAL VARAIBLES:
    integer :: clm_npts
    integer :: clm_surf_npts

    ! Create CLM-PFLOTRAN mapping files
    call CreateCLMPFLOTRANInterfaceDate(this, bounds_clump, clm_npts, clm_surf_npts)

    ! Create CLM-PFLOTRAN mapping files
    call CreateCLMPFLOTRANMaps(this, bounds_clump, clm_npts, clm_surf_npts)

    ! Initialize PFLOTRAN states
    call pflotranModelStepperRunInit(this%pflotran_m)

    ! Get top surface area
    call pflotranModelGetTopFaceArea(this%pflotran_m)

    ! Get PFLOTRAN states
    call pflotranModelGetUpdatedData(this%pflotran_m)

    ! Save the data need by ELM
    call extract_data_for_elm(this, l2e_init_list, e2l_init_list, bounds_clump)

  end subroutine EM_PFLOTRAN_Init

  !-----------------------------------------------------------------------
  subroutine CreateCLMPFLOTRANInterfaceDate(this, bounds, clm_npts, clm_surf_npts)
    !
    ! !DESCRIPTION:
    ! Allocates memory for CLM-PFLOTRAN interface
    !
    ! !USES:
    use decompMod                     , only : bounds_type
    use clm_varpar                    , only : nlevsoi, nlevgrnd
    use shr_kind_mod                  , only: r8 => shr_kind_r8
    ! pflotran
    use Option_module                 , only : printErrMsg
    use Simulation_Base_class         , only : simulation_base_type
    use Simulation_Subsurface_class   , only : simulation_subsurface_type
    use Simulation_Surface_class      , only : simulation_surface_type
    use Simulation_Surf_Subsurf_class , only : simulation_surfsubsurface_type
    use Realization_Subsurface_class  , only : realization_subsurface_type
    use Realization_Surface_class     , only : realization_surface_type
    use PFLOTRAN_Constants_module
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_pflotran_type) , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(inout) :: clm_npts
    integer                 , intent(inout) :: clm_surf_npts
    !
    ! LOCAL VARAIBLES:
    integer                                      :: nlevmapped
    character(len= 128)                          :: subname = 'CreateCLMPFLOTRANInterfaceDate' ! subroutine name
    class(realization_subsurface_type) , pointer :: realization
    class(realization_surface_type)    , pointer :: surf_realization

    ! Initialize PETSc vector for data transfer between CLM and PFLOTRAN
    call CLMPFLOTRANIDataInit()

    nullify(realization)
    nullify(surf_realization)

    select type (simulation => this%pflotran_m%simulation)
    class is (simulation_subsurface_type)
       realization => simulation%realization
       nullify(surf_realization)
    class is (simulation_surface_type)
       nullify(realization)
       surf_realization => simulation%surf_realization
    class is (simulation_surfsubsurface_type)
       realization => simulation%realization
       surf_realization => simulation%surf_realization
    class default
       this%pflotran_m%option%io_buffer = "clm-pflotran only works with surface and subsurface simulations."
       write(*, '(/A/)') this%pflotran_m%option%io_buffer
       call printErrMsg(this%pflotran_m%option)
    end select

    ! Compute number of cells in CLM domain.
    ! Assumption-1: One column per CLM grid cell.

    ! Check if the number of CLM vertical soil layers defined in the mapping
    ! file read by PFLOTRAN matches either nlevsoi or nlevgrnd
    clm_pf_idata%nzclm_mapped = this%pflotran_m%map_clm_sub_to_pf_sub%clm_nlevsoi
    nlevmapped                = clm_pf_idata%nzclm_mapped
    if ( (nlevmapped /= nlevsoi) .and. (nlevmapped /= nlevgrnd) ) then
       call endrun(trim(subname)//' ERROR: Number of layers PFLOTRAN thinks CLM should '// &
            'have do not match either nlevsoi or nlevgrnd. Abortting' )
    end if

    clm_npts = (bounds%endg - bounds%begg + 1)*nlevmapped
    clm_surf_npts = (bounds%endg - bounds%begg + 1)

    ! CLM: Subsurface domain (local and ghosted cells)
    clm_pf_idata%nlclm_sub = clm_npts
    clm_pf_idata%ngclm_sub = clm_npts

    ! CLM: Surface of subsurface domain (local and ghosted cells)
    clm_pf_idata%nlclm_2dsub = (bounds%endg - bounds%begg + 1)
    clm_pf_idata%ngclm_2dsub = (bounds%endg - bounds%begg + 1)
    ! For CLM: Same as surface of subsurface domain
    clm_pf_idata%nlclm_srf = clm_surf_npts
    clm_pf_idata%ngclm_srf = clm_surf_npts

    ! PFLOTRAN: Subsurface domain (local and ghosted cells)
    clm_pf_idata%nlpf_sub = realization%patch%grid%nlmax
    clm_pf_idata%ngpf_sub = realization%patch%grid%ngmax

    ! PFLOTRAN: Surface of subsurface domain (local and ghosted cells)
    if (this%pflotran_m%option%iflowmode == TH_MODE) then
       clm_pf_idata%nlpf_2dsub = pflotranModelNSurfCells3DDomain(this%pflotran_m)
       clm_pf_idata%ngpf_2dsub = pflotranModelNSurfCells3DDomain(this%pflotran_m)
    else
       clm_pf_idata%nlpf_2dsub = 0
       clm_pf_idata%ngpf_2dsub = 0
    endif

    ! PFLOTRAN: Surface domain (local and ghosted cells)
    if (associated(surf_realization) .and. this%pflotran_m%option%nsurfflowdof > 0) then
       clm_pf_idata%nlpf_srf = surf_realization%patch%grid%nlmax
       clm_pf_idata%ngpf_srf = surf_realization%patch%grid%ngmax
    else
       clm_pf_idata%nlpf_srf = 0
       clm_pf_idata%ngpf_srf = 0
    endif

    ! Allocate vectors for data transfer between CLM and PFLOTRAN.
    call CLMPFLOTRANIDataCreateVec(MPI_COMM_WORLD)

  end subroutine CreateCLMPFLOTRANInterfaceDate

  !-----------------------------------------------------------------------
  subroutine CreateCLMPFLOTRANMaps(this, bounds, clm_npts, clm_surf_npts)
    !
    ! !DESCRIPTION:
    ! Creates maps to transfer data between CLM and PFLOTRAN
    !
    ! !USES:
    use clm_varctl      , only : pflotran_surfaceflow, pflotran_th_mode, pflotran_th_freezing
    use decompMod       , only : bounds_type, ldecomp
    use PFLOTRAN_Constants_module
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_pflotran_type) , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(inout) :: clm_npts
    integer                 , intent(inout) :: clm_surf_npts
    !
    ! LOCAL VARAIBLES:
    integer             :: g,j
    integer             :: nlevmapped
    integer, pointer    :: clm_cell_ids_nindex(:)
    integer, pointer    :: clm_surf_cell_ids_nindex(:)

    ! Save cell IDs of CLM grid
    allocate(clm_cell_ids_nindex(     1:clm_npts     ))
    allocate(clm_surf_cell_ids_nindex(1:clm_surf_npts))

    nlevmapped     = clm_pf_idata%nzclm_mapped
    clm_npts       = 0
    clm_surf_npts  = 0
    do g = bounds%begg, bounds%endg
       do j = 1,nlevmapped
          clm_npts = clm_npts + 1
          clm_cell_ids_nindex(clm_npts) = (ldecomp%gdc2glo(g)-1)*nlevmapped + j - 1
       enddo
       clm_surf_npts = clm_surf_npts + 1
       clm_surf_cell_ids_nindex(clm_surf_npts) = (ldecomp%gdc2glo(g)-1)*nlevmapped
    enddo

    ! Initialize maps for transferring data between CLM and PFLOTRAN.
    call pflotranModelInitMapping(this%pflotran_m, clm_cell_ids_nindex, &
                                  clm_npts, CLM_SUB_TO_PF_SUB)
    call pflotranModelInitMapping(this%pflotran_m, clm_cell_ids_nindex, &
                                  clm_npts, CLM_SUB_TO_PF_EXTENDED_SUB)
    call pflotranModelInitMapping(this%pflotran_m, clm_cell_ids_nindex, &
                                  clm_npts, PF_SUB_TO_CLM_SUB)

    if (this%pflotran_m%option%iflowmode == TH_MODE) then
      pflotran_th_mode = .true.
      if (this%pflotran_m%option%use_th_freezing) pflotran_th_freezing  = .true.
    endif

    if (this%pflotran_m%option%nsurfflowdof > 0) then
      pflotran_surfaceflow = .true.
      call pflotranModelInitMapping(this%pflotran_m, clm_surf_cell_ids_nindex, &
                                    clm_surf_npts, PF_SRF_TO_CLM_SRF)
      call pflotranModelInitMapping(this%pflotran_m, clm_surf_cell_ids_nindex, &
                                    clm_surf_npts, CLM_SRF_TO_PF_SRF)
    else
      if (this%pflotran_m%option%iflowmode == TH_MODE) then
        call pflotranModelInitMapping(this%pflotran_m, clm_surf_cell_ids_nindex, &
                                      clm_surf_npts, CLM_SRF_TO_PF_2DSUB)
      endif
    endif

    deallocate(clm_cell_ids_nindex)
    deallocate(clm_surf_cell_ids_nindex)

  end subroutine CreateCLMPFLOTRANMaps


  !-----------------------------------------------------------------------
  subroutine extract_data_for_elm(this, l2e_init_list, e2l_init_list, bounds_clump)
    !
    !DESCRIPTION
    !  Saves
    !
    use landunit_varcon           , only : istsoil
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_MASS
    use MultiPhysicsProbConstants , only : VAR_SOIL_MATRIX_POT
    use clm_varcon                , only : denice, denh2o
    use clm_varpar                , only : nlevgrnd
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_pflotran_type)              :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    type(bounds_type)    , intent(in)    :: bounds_clump

    ! !LOCAL VARIABLES:
    integer               :: c,j,g,l,pf_j  ! do loop indices
    integer               :: nlevmapped
    integer               :: gcount
    integer               :: bounds_proc_begc, bounds_proc_endc

    integer     , pointer :: col_active(:)
    integer     , pointer :: col_landunit(:)
    integer     , pointer :: col_gridcell(:)
    integer     , pointer :: lun_type(:)

    real(r8)    , pointer :: l2e_h2osoi_ice(:,:)

    real(r8)    , pointer :: e2l_h2osoi_liq(:,:)
    real(r8)    , pointer :: e2l_h2osoi_ice(:,:)
    real(r8)    , pointer :: e2l_h2osoi_vol(:,:)
    real(r8)    , pointer :: e2l_zwt(:)
    real(r8)    , pointer :: e2l_mflx_snowlyr_col(:)
    real(r8)    , pointer :: e2l_watsatc(:,:)
    real(r8)    , pointer :: e2l_hksatc(:,:)
    real(r8)    , pointer :: e2l_bswc(:,:)
    real(r8)    , pointer :: e2l_sucsatc(:,:)

    real(r8)    , pointer :: dz(:,:)

    PetscScalar , pointer :: sat_clm_loc(:)
    PetscScalar , pointer :: watsat_clm_loc(:)
    PetscScalar , pointer :: hksat_clm_loc(:)
    PetscScalar , pointer :: bsw_clm_loc(:)
    PetscScalar , pointer :: sucsat_clm_loc(:)
    PetscErrorCode        :: ierr

    character(len= 128)   :: subname = 'extract_data_for_elm'
    !-----------------------------------------------------------------------

    bounds_proc_begc = bounds_clump%begc
    bounds_proc_endc = bounds_clump%endc
    nlevmapped       = clm_pf_idata%nzclm_mapped

    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_active             , col_active   )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_landunit_index     , col_landunit )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_gridcell_index     , col_gridcell )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_dz                , dz           )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_landunit_type          , lun_type     )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_h2osoi_ice            , l2e_h2osoi_ice       )

    call e2l_init_list%GetPointerToReal1D(this%index_e2l_init_state_wtd             , e2l_zwt              )
    call e2l_init_list%GetPointerToReal1D(this%index_e2l_init_flux_mflx_snowlyr_col , e2l_mflx_snowlyr_col )

    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_liq      , e2l_h2osoi_liq       )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_ice      , e2l_h2osoi_ice       )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_vol      , e2l_h2osoi_vol       )

    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_watsatc     , e2l_watsatc )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_hksatc      , e2l_hksatc  )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_bswc        , e2l_bswc    )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_sucsatc     , e2l_sucsatc )
    
    call pflotranModelGetSoilProp(this%pflotran_m)

    ! Set initial value of for ELM
    e2l_mflx_snowlyr_col(:) = 0._r8
    e2l_zwt(:)              = 0._r8

    ! Initialize soil moisture
    call VecGetArrayF90(clm_pf_idata%sat_clm      , sat_clm_loc    , ierr)
    call VecGetArrayF90(clm_pf_idata%watsat2_clm  , watsat_clm_loc , ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_x2_clm , hksat_clm_loc  , ierr)
    call VecGetArrayF90(clm_pf_idata%bsw2_clm     , bsw_clm_loc    , ierr)
    call VecGetArrayF90(clm_pf_idata%sucsat2_clm  , sucsat_clm_loc , ierr)

    do c = bounds_proc_begc, bounds_proc_endc
       if (col_active(c) == 1) then
          l = col_landunit(c)
          if (lun_type(l) == istsoil) then
             g = col_gridcell(c)
             gcount = g - bounds_clump%begg
             do j = 1, nlevgrnd
                pf_j = gcount*nlevmapped + j

                if (j <= nlevmapped) then
                   e2l_h2osoi_liq(c,j) = sat_clm_loc(pf_j)*watsat_clm_loc(pf_j)*dz(c,j)*1.e3_r8

                   e2l_h2osoi_vol(c,j) = e2l_h2osoi_liq(c,j)/dz(c,j)/denh2o + &
                        l2e_h2osoi_ice(c,j)/dz(c,j)/denice
                   e2l_h2osoi_vol(c,j) = min(e2l_h2osoi_vol(c,j),watsat_clm_loc(pf_j))
                   e2l_h2osoi_ice(c,j) = 0._r8

                   e2l_watsatc(c,j) = watsat_clm_loc(pf_j)
                   e2l_hksatc(c,j)  = hksat_clm_loc(pf_j)
                   e2l_bswc(c,j)    = bsw_clm_loc(pf_j)
                   e2l_sucsatc(c,j) = sucsat_clm_loc(pf_j)
                else
                   e2l_h2osoi_liq(c,j) = e2l_h2osoi_liq(c,nlevmapped)
                   e2l_h2osoi_vol(c,j) = e2l_h2osoi_vol(c,nlevmapped)
                   e2l_h2osoi_ice(c,j) = 0._r8
                   e2l_watsatc(c,j)    = e2l_watsatc(c,nlevmapped)
                   e2l_hksatc(c,j)     = e2l_hksatc(c,nlevmapped)
                   e2l_bswc(c,j)       = e2l_bswc(c,nlevmapped)
                   e2l_sucsatc(c,j)    = e2l_sucsatc(c,nlevmapped)
                end if

             enddo
          else
             write(iulog,*)'WARNING: Land Unit type other than soil type is present within the domain'
             call endrun( trim(subname)//' ERROR: Land Unit type not supported' )             
          endif
       endif
    enddo

    call VecRestoreArrayF90(clm_pf_idata%sat_clm      , sat_clm_loc    , ierr)
    call VecRestoreArrayF90(clm_pf_idata%watsat2_clm  , watsat_clm_loc , ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_x2_clm , hksat_clm_loc  , ierr)
    call VecRestoreArrayF90(clm_pf_idata%bsw2_clm     , bsw_clm_loc    , ierr)
    call VecRestoreArrayF90(clm_pf_idata%sucsat2_clm  , sucsat_clm_loc , ierr)

   end subroutine extract_data_for_elm

   !------------------------------------------------------------------------
   subroutine EM_PFLOTRAN_Solve(this, em_stage, dt, nstep, clump_rank, l2e_list, e2l_list, &
        bounds_clump)
    !
    ! !DESCRIPTION:
    ! The VSFM dirver subroutine
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_pflotran_type)              :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump


  end subroutine EM_PFLOTRAN_Solve

#endif
end module ExternalModelPFLOTRANMod
