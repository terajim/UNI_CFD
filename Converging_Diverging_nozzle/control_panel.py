import conservative_sub_super
import non_conservative_sub_super
import subsonic_conservative
import subsonic_non_conservative
import sock_capture_conservative
import sock_capture_non_conservative

### Αυτό είναι το κεντρικό control panel. Ενεργοποιόντας οποιαδήποτε από τις εξισώσεις,
### θα τρέξει και η αντίστοιχη επίλυσης της ροής
### Νοte: Το sock capture non conservative είναι προσωρινά oyt of order



#conservative_sub_super.conservative(80, 0.75, 5000, 10**(-5))

#non_conservative_sub_super.nonconservative(41, 0.75, 5000, 10**(-5))

#subsonic_non_conservative.subsonic(0.6, 0.93, 41, 10**(-4), 10000)

#subsonic_conservative.conservative(41, 0.6, 10000, 0.93, 10**(-4))

sock_capture_conservative.conservative(41, 0.5, 10000, 0.6784, 0.3, 10**(-4))

#sock_capture_non_conservative.noncons(0.5, 0.6784, 41, 10**(-4),10000, 0.1, 0.1)
