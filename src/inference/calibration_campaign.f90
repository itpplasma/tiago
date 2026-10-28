module tiago_calibration_campaign
    use fortnum_gaussian_constraints, only: campaign_factor => brownian_factor, &
                                            constrain_discrepancy
    implicit none
    private
    public :: campaign_factor, constrain_discrepancy
end module
