%%  "setvar" sets variables or options as global variable
%%
%%  [input]
%%   (stype, ...)
%%
%%   ('data', my, nlag) -> set data and lags
%%
%%   ('intercept', fli) -> set time-varying intercept
%%      0: off, demean (default)
%%      1: on
%%
%%  ('LTb', flLT) -> set threshold for b
%%      0: off
%%      1: on (default)
%%  ('LTa', flLT) -> set threshold for a
%%      0: off
%%      1: on (default)
%%
%%   ('ranseed', iseed) -> set ranseed (default: 1)
%%
%%  [global variables]
%%      m_my:     variable y_t (n*k matrix)
%%      m_asvar:  variable names
%%      m_nl:     # of lags (l, scalar)
%%      m_ns:     # of time periods (n, scalar)
%%      m_nk:     # of series (k, scalar)
%%      m_fli:    intercept type (flag)
%%      m_flLTb:  latent threshold for b (flag)
%%      m_flLTa:  latent threshold for a (flag)
%%      m_iseed:  ranseed


function lt_setvar(stype, arg1, arg2, arg3)

global m_my m_asvar m_nl m_ns m_nk m_fli m_flLTb m_flLTa m_iseed

switch stype
    case 'data'
        m_my = arg1;
        m_asvar = arg2;
        m_nl = arg3;
        m_ns = size(m_my, 1);
        m_nk = size(m_my, 2);
    case 'intercept'
        m_fli = arg1;
    case 'LTb'
        m_flLTb = arg1;
    case 'LTa'
        m_flLTa = arg1;
    case 'ranseed'
        m_iseed = arg1;
end
