/*
    //this is a very un-scala function.... + no exception handling
    def modpLgndr(l: Int, m: Int, x: Double): Double = {
      assert(0 <= m && m <= l && abs(x) <= 1.0)
      val dl = l.toDouble
      val dm = m.toDouble
      val norm = sqrt(2.0 * dl + 1.0) / sqrt(4.0 * Pi)
      var pmm = norm
      if (m != 0) pmm = (pow(-1, m)).toDouble * pmm * xfact(2 * m) * pow((1.0-x * x), (dm / 2.0))
      if (l == m) pmm
      else {
        var pmmp1 = x * pmm * sqrt(2.0 * m + 1.0)
        if (l == m + 1) pmmp1
        else {
          var pll = 0.0
          var dll = 0.0
          for (ll <- m + 2 to l) {
            dll = ll.toDouble
            pll = (x * (2.0 * dll - 1.0) * pmmp1 - sqrt(pow((dll - 1.0), 2.0) - dm * dm) * pmm) / sqrt(pow(dll, 2.0) - pow(dm, 2.0))
            pmm = pmmp1
            pmmp1 = pll
          }
          pll
        }
      }
    }
    */
