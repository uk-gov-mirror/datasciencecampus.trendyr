##' Outputs a connected scatterplot
##'
##' @param df The data frame with date, trend and first derivative
##' @param label_type The date format of the time series. There are four options: ymd, md, my, y.
##'
##' @author Sonia Mazzi and Harry Churchley
##' @export
flightpath <- function(df = randomData, label_type = labeltype) {
  aux <- df %>%
    mutate(year = as.character(year(date)), month = month(date, label=T), day = as.character(day(date))) %>%
    mutate(label = dplyr::case_when(label_type == "ymd" ~ paste0(month, " ", day, ", ", year),
                                    label_type == "md" ~ paste0(month, " ", day),
                                    label_type == "my" ~ paste0(month, " ", year),
                                    label_type == "y" ~ year))

  pp <- aux %>%
    ggplot(aes(x = trend, y = first_derivative, label = label)) +
    geom_point() +
    geom_text_repel(data = sample_frac(aux, 0.3), size = 2.8) +
    geom_segment(color = "#69b3a2",
                 aes(xend = c(tail(trend, n = -1), NA),
                     yend = c(tail(first_derivative, n = -1), NA)
                 ),
                 arrow = arrow(length = unit(0.3, "cm"))
    ) +
    #ggtitle("title") +
    geom_hline(yintercept = 0, colour = "orange") +
    labs(y = "First derivative", x = "Trend") +

    scale_colour_manual(name = " ", values = cols) +
    theme(legend.position = "bottom") +
    theme(legend.key = element_blank())

  return(pp)
}
